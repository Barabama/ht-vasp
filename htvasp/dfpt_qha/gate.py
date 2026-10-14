"""The VASP input gate: check the VALUE of every tag, not the presence of a key.

Why this exists
---------------
Three consecutive incidents in this project had the same shape: a builder wrote
a *template* and the check asked only whether *a key was present*.

  1. POSCAR lattice rows written as literal ``"{lat:.16f}"`` strings -> VASP read
     "0 types and 0 ions" and then failed with an unrelated LUSE_VDW message.
  2. INCAR got ``ISMEAR = ISMEAR`` because a dict was unpacked as a 2-tuple, and
     20 directories "passed" a check that only looked for the substring.
  3. A 16-atom cell got the 2-atom MAGMOM -> VASP refused with "Error reading
     item MAGMOM from file INCAR. Error code was IERR= 5", killing 14 cases.
     The validator had checked POTCAR order, atom counts, POSCAR numerics,
     KPOINTS and ISMEAR -- everything except MAGMOM length.

So the contract here is: check the value of every tag, **against the structure
the directory actually contains**.  Every problem is reported as what was read
and what the contract requires.

Two layers
----------
``validate_dir``   static: pure text and structure checks, no VASP needed.
``dryrun_case``    dynamic: copy to scratch, run one electronic step, and
                   confirm VASP itself echoed NIONS / NKPTS / ISPIN without an
                   "Error reading item".  This is the only check that models
                   VASP's own parser, which is where the MAGMOM-length and
                   POSCAR-layout failures came from.  It needs a compute node,
                   so it is invoked from a submitted job.

The gate returns a ``GateResult`` rather than mutating a dict.  Callers that
need to attach their own metadata to a ``case.json`` compose it explicitly::

    gate = validate_dir(case_dir)
    record = {**gate.model_dump(mode="json"), "kgrid": n, **meta}

which keeps the verdict itself immutable -- a gate whose result can be edited
after the fact is not a gate.
"""

from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Sequence

from pydantic import BaseModel, ConfigDict, Field

from htvasp.dfpt_qha.recipe import axis_contains_gamma, gamma_centred

# Tags whose value must be an integer / a float.  A tag present with a
# non-literal value (e.g. "ISMEAR = ISMEAR") is a failure, not a warning.
INT_TAGS = (
    "ISMEAR", "ISPIN", "IBRION", "NSW", "ISIF", "NELM", "NELMIN", "NCORE",
    "KPAR", "NFREE", "ISYM", "ISTART", "ICHARG", "LMAXMIX", "NELMDL",
)
FLOAT_TAGS = ("SIGMA", "ENCUT", "EDIFF", "EDIFFG", "POTIM", "SYMPREC")
BOOL_TAGS = ("LWAVE", "LCHARG", "LREAL", "LOPTICS", "LVTOT", "LDAU", "LASPH")

INCAR_NAMES = ("INCAR",)

# VASP's own failure strings.  Matching these is how the dry run learns that
# VASP rejected an input rather than merely producing no OUTCAR.
VASP_REJECTIONS = (
    "Error reading item",
    "I REFUSE TO CONTINUE",
    "POSCAR found : 0 types",
    "not enough POTCAR",
)

_SHIFT_TOL = 1e-9


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()[:16]


_REPEAT = re.compile(r"^(\d+)\*(.+)$")


def _same_element(potcar_element: str, poscar_species: str) -> bool:
    """Whether one POTCAR block is the potential for one POSCAR species.

    Not a string equality.  potpaw_PBE has no plain ``Nb`` -- only ``Nb_pv``
    and ``Nb_sv`` -- so the pool's Nb endpoints legitimately carry a POTCAR
    whose TITEL reads ``PAW_PBE Nb_pv`` against a POSCAR species of ``Nb``.  A
    strict compare rejects every one of them.

    The suffix must follow an underscore, so ``N`` still does not match ``Nb``.
    """

    potcar_element, poscar_species = potcar_element.upper(), poscar_species.upper()
    return potcar_element == poscar_species or potcar_element.startswith(
        f"{poscar_species}_"
    )


def parse_repeated_floats(text: str) -> list[float] | None:
    """Values from an INCAR list tag, expanding VASP's ``n*value`` syntax.

    VASP accepts ``MAGMOM = 2*0.0`` for two zero moments, and pymatgen writes
    exactly that whenever the values repeat.  A reader that treats the tag as a
    plain token list gets this wrong twice over: it reports ``'2*0.0'`` as a
    non-numeric value, and it counts a run of N as ONE entry -- so a MAGMOM of
    the wrong length can look like the right length.

    Returns None if any value is not a number.
    """

    values: list[float] = []
    for token in text.split():
        match = _REPEAT.match(token)
        if match is None:
            try:
                values.append(float(token))
            except ValueError:
                return None
            continue
        count = int(match.group(1))
        try:
            value = float(match.group(2))
        except ValueError:
            return None
        if count < 1:
            return None
        values.extend([value] * count)
    return values


class PoscarInfo(BaseModel):
    """What the POSCAR actually says.  ``error`` is set instead of the rest."""

    model_config = ConfigDict(extra="forbid")

    error: str | None = None
    lines: int = 0
    scale: float | None = None
    lattice: list[list[float]] = Field(default_factory=list)
    species: list[str] = Field(default_factory=list)
    counts: list[int] = Field(default_factory=list)
    natoms: int = 0
    coord_mode: str | None = None
    n_coord_lines: int = 0
    #: species in POSCAR order, deduplicated -- the order VASP matches POTCAR by
    species_order: list[str] = Field(default_factory=list)
    #: one entry per atom, so a MAGMOM expansion can be checked against it
    per_atom_species: list[str] = Field(default_factory=list)
    volume_angstrom3: float = 0.0
    sha256: str | None = None


class PotcarInfo(BaseModel):
    model_config = ConfigDict(extra="forbid")

    error: str | None = None
    blocks: int = 0
    titles: list[str] = Field(default_factory=list)
    elements: list[str] = Field(default_factory=list)
    zval: list[float] = Field(default_factory=list)
    enmax: list[float] = Field(default_factory=list)
    sha256: str | None = None


class KpointsInfo(BaseModel):
    model_config = ConfigDict(extra="forbid")

    error: str | None = None
    mode: str | None = None
    scheme_line: str | None = None
    mesh: list[int] | None = None
    shift: list[float] | None = None
    #: None means the shift is not one the rule covers, so membership is
    #: UNKNOWN -- not False.  Reporting "excludes Gamma" for an unmodelled
    #: shift would be a guess dressed as a check.
    gamma_centered: bool | None = None
    gamma_rule: str | None = None
    sha256: str | None = None


class IncarInfo(BaseModel):
    model_config = ConfigDict(extra="forbid")

    error: str | None = None
    tags: dict[str, str] = Field(default_factory=dict)
    sha256: str | None = None


class DryRunResult(BaseModel):
    """What VASP itself said about the inputs, after one electronic step."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    dir: str
    ok: bool
    problems: tuple[str, ...] = ()
    scratch: str | None = None
    nions: int | None = None
    nkpts: int | None = None
    ispin: int | None = None
    magmom_echo: int | None = None


class GateResult(BaseModel):
    """The verdict for one VASP directory, plus the reads it was based on.

    Frozen: a gate whose result can be edited after the fact is not a gate, and
    the builders that consume this attach their own metadata to a copy rather
    than mutating the verdict.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    dir: str
    ok: bool
    problems: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    natoms: int | None = None
    species_order: tuple[str, ...] = ()
    volume_angstrom3: float | None = None
    poscar_sha256: str | None = None
    potcar_sha256: str | None = None
    incar_sha256: str | None = None
    kpoints_sha256: str | None = None
    potcar_zval: tuple[float, ...] = ()
    mesh: tuple[int, ...] | None = None
    shift: tuple[float, ...] | None = None
    gamma_centered: bool | None = None
    ismear: str | None = None
    sigma: str | None = None
    #: only populated by the dry run; the static gate never sets it
    dryrun: DryRunResult | None = None


# --------------------------------------------------------------------------
# readers -- each returns what it found, or an ``error``
# --------------------------------------------------------------------------


def read_poscar(directory: Path) -> PoscarInfo:
    """Parse just enough POSCAR to know atom counts, species and order."""

    path = directory / "POSCAR"
    if not path.is_file():
        return PoscarInfo(error="no POSCAR")
    lines = path.read_text().splitlines()
    if len(lines) < 8:
        return PoscarInfo(error=f"POSCAR has only {len(lines)} lines")
    try:
        scale = float(lines[1].split()[0])
    except (ValueError, IndexError):
        return PoscarInfo(error=f"POSCAR line 2 is not a scale factor: {lines[1]!r}")
    lattice: list[list[float]] = []
    for i in (2, 3, 4):
        try:
            row = [float(x) for x in lines[i].split()]
        except ValueError:
            return PoscarInfo(error=f"POSCAR line {i + 1} is not numeric: {lines[i]!r}")
        if len(row) != 3:
            return PoscarInfo(error=f"POSCAR line {i + 1} has {len(row)} numbers, not 3")
        lattice.append(row)
    index = 5
    if lines[index].strip()[:1].lower() == "s":  # optional "Selective dynamics"
        index += 1
    species = lines[index].split()
    index += 1
    try:
        counts = [int(x) for x in lines[index].split()]
    except ValueError:
        return PoscarInfo(error=f"POSCAR count line is not integers: {lines[index]!r}")
    if len(counts) != len(species):
        return PoscarInfo(
            error=f"POSCAR lists {len(species)} species but {len(counts)} counts"
        )
    mode = lines[index + 1].strip()[:1].lower()
    body = [line for line in lines[index + 2:] if line.strip()]
    a, b, c = lattice
    volume = (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    return PoscarInfo(
        lines=len(lines),
        scale=scale,
        lattice=lattice,
        species=species,
        counts=counts,
        natoms=sum(counts),
        coord_mode="Direct" if mode == "d" else ("Cartesian" if mode == "c" else mode),
        n_coord_lines=len(body),
        species_order=list(dict.fromkeys(species)),
        per_atom_species=[s for s, c in zip(species, counts) for _ in range(c)],
        volume_angstrom3=abs(volume) * abs(scale) ** 3,
        sha256=_sha(path.read_bytes()),
    )


def read_potcar(directory: Path) -> PotcarInfo:
    path = directory / "POTCAR"
    if not path.is_file():
        return PotcarInfo(error="no POTCAR")
    text = path.read_text(errors="replace")
    titles = re.findall(r"TITEL\s*=\s*(\S+)\s+(\w+)", text)
    return PotcarInfo(
        blocks=len(titles),
        titles=[t[0] for t in titles],
        elements=[t[1] for t in titles],
        zval=[float(z) for z in re.findall(r"ZVAL\s*=\s*([\d.]+)", text)],
        enmax=[float(e) for e in re.findall(r"ENMAX\s*=\s*([\d.]+)", text)],
        sha256=_sha(path.read_bytes()),
    )


def read_kpoints(directory: Path) -> KpointsInfo:
    path = directory / "KPOINTS"
    if not path.is_file():
        return KpointsInfo(error="no KPOINTS")
    lines = [line for line in path.read_text().splitlines() if line.strip() != ""]
    if len(lines) < 3:
        return KpointsInfo(error=f"KPOINTS has only {len(lines)} non-blank lines")
    scheme = lines[2].strip()
    digest = _sha(path.read_bytes())
    if "monkhorst" not in scheme.lower() and "gamma" not in scheme.lower():
        return KpointsInfo(
            mode="explicit",
            scheme_line=scheme,
            gamma_rule="explicit list: Gamma membership must be checked in IBZKPT",
            sha256=digest,
        )
    try:
        mesh = [int(x) for x in lines[3].split()]
    except (ValueError, IndexError):
        return KpointsInfo(error=f"KPOINTS mesh line is not integers: {lines[3]!r}")
    if len(mesh) != 3:
        return KpointsInfo(error=f"KPOINTS mesh has {len(mesh)} entries, not 3")
    if len(lines) >= 5:
        try:
            shift = [float(x) for x in lines[4].split()]
        except ValueError:
            return KpointsInfo(error=f"KPOINTS shift line is not numeric: {lines[4]!r}")
    else:
        shift = [0.0, 0.0, 0.0]
    if len(shift) != 3:
        return KpointsInfo(error=f"KPOINTS shift has {len(shift)} entries, not 3")

    uniform_zero = all(abs(x) < _SHIFT_TOL for x in shift)
    uniform_half = all(abs(x - 0.5) < _SHIFT_TOL for x in shift)
    if uniform_zero or uniform_half:
        label = "shift 0" if uniform_zero else "shift 0.5"
        centred = gamma_centred(tuple(mesh), tuple(shift))
        offenders = [
            f"axis {i} (n={n})"
            for i, (n, s) in enumerate(zip(mesh, shift))
            if not axis_contains_gamma(n, s)
        ]
        rule = f"{label} with mesh {mesh} " + (
            "includes Gamma"
            if centred
            else f"excludes Gamma at {', '.join(offenders)}; "
            "the folding-compatible family is 'even n + shift 0.5'"
        )
    else:
        # Membership is genuinely not modelled for a mixed or arbitrary shift,
        # so it is reported as unknown rather than as a failure.
        centred, rule = None, f"non-standard shift {shift}: Gamma membership not modelled"
    return KpointsInfo(
        mode="automatic",
        scheme_line=scheme,
        mesh=mesh,
        shift=shift,
        gamma_centered=centred,
        gamma_rule=rule,
        sha256=digest,
    )


def read_incar(directory: Path) -> IncarInfo:
    for name in INCAR_NAMES:
        path = directory / name
        if path.is_file():
            break
    else:
        return IncarInfo(error="no INCAR")
    tags: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.split("#")[0].split("!")[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        tags[key.strip().upper()] = value.strip()
    return IncarInfo(tags=tags, sha256=_sha(path.read_bytes()))


# --------------------------------------------------------------------------
# the static contract
# --------------------------------------------------------------------------


def validate_dir(
    directory: Path,
    *,
    expect_natoms: int | None = None,
    expect_species: Sequence[str] | None = None,
) -> GateResult:
    """Check one VASP directory.  Returns the verdict; never raises on content."""

    problems: list[str] = []
    warnings: list[str] = []

    poscar = read_poscar(directory)
    potcar = read_potcar(directory)
    kpoints = read_kpoints(directory)
    incar = read_incar(directory)

    for name, read in (
        ("POSCAR", poscar), ("POTCAR", potcar),
        ("KPOINTS", kpoints), ("INCAR", incar),
    ):
        if read.error:
            problems.append(f"{name}: {read.error}")
    if problems:
        return GateResult(dir=str(directory), ok=False, problems=problems)

    natoms = poscar.natoms

    # ---- POSCAR internals -------------------------------------------------
    if poscar.n_coord_lines != natoms:
        problems.append(
            f"POSCAR declares {natoms} atoms but lists "
            f"{poscar.n_coord_lines} coordinate lines"
        )
    if expect_natoms is not None and natoms != expect_natoms:
        problems.append(f"POSCAR has {natoms} atoms, expected {expect_natoms}")
    if expect_species is not None and [s.upper() for s in poscar.species_order] != [
        s.upper() for s in expect_species
    ]:
        problems.append(
            f"POSCAR species order {poscar.species_order} != expected "
            f"{list(expect_species)}"
        )
    if any(count == 0 for count in poscar.counts):
        problems.append(f"POSCAR has a zero species count: {poscar.counts}")
    if abs(poscar.scale or 0.0) == 0:
        problems.append("POSCAR scale factor is 0")
    if poscar.volume_angstrom3 <= 0:
        problems.append(f"POSCAR volume is not positive: {poscar.volume_angstrom3}")

    # ---- POTCAR vs POSCAR -------------------------------------------------
    # VASP matches POTCAR to POSCAR POSITIONALLY, so the order is the contract.
    if potcar.blocks != len(poscar.species):
        problems.append(
            f"POTCAR has {potcar.blocks} blocks but POSCAR declares "
            f"{len(poscar.species)} species"
        )
    if len(potcar.elements) == len(poscar.species) and not all(
        _same_element(p, s) for p, s in zip(potcar.elements, poscar.species)
    ):
        problems.append(
            f"POTCAR element order {potcar.elements} != POSCAR order "
            f"{poscar.species} (VASP matches them positionally)"
        )
    if len(potcar.zval) != potcar.blocks:
        problems.append(
            f"POTCAR has {potcar.blocks} blocks but {len(potcar.zval)} ZVAL values"
        )

    # ---- KPOINTS vs POSCAR ------------------------------------------------
    if kpoints.mode == "automatic":
        if min(kpoints.mesh or [0]) < 1:
            problems.append(f"KPOINTS mesh {kpoints.mesh} has a non-positive entry")
        if kpoints.gamma_centered is None:
            warnings.append(f"KPOINTS {kpoints.gamma_rule}")
        elif not kpoints.gamma_centered:
            warnings.append(f"KPOINTS {kpoints.gamma_rule}")
    else:
        warnings.append(
            "KPOINTS is an explicit list; Gamma membership must be read from IBZKPT"
        )

    # ---- INCAR value contracts -------------------------------------------
    tags = incar.tags
    for tag in INT_TAGS:
        if tag in tags:
            try:
                int(tags[tag])
            except ValueError:
                problems.append(f"INCAR {tag} = {tags[tag]!r} is not an integer literal")
    for tag in FLOAT_TAGS:
        if tag in tags:
            try:
                float(tags[tag])
            except ValueError:
                problems.append(f"INCAR {tag} = {tags[tag]!r} is not numeric")
    for tag in BOOL_TAGS:
        if tag in tags and tags[tag].strip(".").upper() not in ("TRUE", "FALSE", "T", "F"):
            problems.append(f"INCAR {tag} = {tags[tag]!r} is not a logical value")

    def as_int(tag: str, default: int = 0) -> int | None:
        """int value of a tag, or None if absent or unparseable.

        The INT_TAGS loop already reports unparseable values; this must not
        raise a second time on the same bad input.
        """

        if tag not in tags:
            return default
        try:
            return int(tags[tag])
        except ValueError:
            return None

    ispin = as_int("ISPIN", 1)
    if "MAGMOM" in tags:
        values = parse_repeated_floats(tags["MAGMOM"])
        if values is None:
            problems.append(f"INCAR MAGMOM has a non-numeric value: {tags['MAGMOM']!r}")
            values = []
        # VASP accepts only one value or exactly NIONS; this is the check whose
        # absence killed 14 cases.
        if values and len(values) not in (1, natoms):
            problems.append(
                f"INCAR MAGMOM has {len(values)} values but the cell has {natoms} "
                f"atoms (VASP accepts only 1 or NIONS)"
            )
        if ispin == 1:
            warnings.append("ISPIN=1 with MAGMOM set: MAGMOM is ignored")
    elif ispin == 2:
        warnings.append("ISPIN=2 with no MAGMOM: VASP defaults to 1.0 per atom")

    if "NCORE" in tags and "KPAR" in tags:
        problems.append(
            "INCAR sets both NCORE and KPAR; VASP documents them as mutually "
            "exclusive parallelisation controls"
        )

    icharg = as_int("ICHARG", 0)
    present = {f.name for f in directory.iterdir()}
    if icharg in (1, 11) and "CHGCAR" not in present:
        problems.append(f"ICHARG={tags['ICHARG']} needs a CHGCAR in the directory; none present")
    if icharg == 1 and as_int("ISTART", 0) != 1:
        warnings.append(f"ICHARG=1 with ISTART={tags.get('ISTART', '0')}; VASP resets ISTART to 1")

    ib_rion = as_int("IBRION", -1)
    if ib_rion is not None and ib_rion not in (-1, 0, 1, 2, 3, 5, 6, 7, 8):
        warnings.append(f"IBRION={ib_rion} is unusual")

    return GateResult(
        dir=str(directory),
        ok=not problems,
        problems=problems,
        warnings=warnings,
        natoms=natoms,
        species_order=poscar.species_order,
        volume_angstrom3=poscar.volume_angstrom3,
        poscar_sha256=poscar.sha256,
        potcar_sha256=potcar.sha256,
        incar_sha256=incar.sha256,
        kpoints_sha256=kpoints.sha256,
        potcar_zval=potcar.zval,
        mesh=kpoints.mesh,
        shift=kpoints.shift,
        gamma_centered=kpoints.gamma_centered,
        ismear=tags.get("ISMEAR"),
        sigma=tags.get("SIGMA"),
    )


# --------------------------------------------------------------------------
# the dynamic gate: let VASP itself parse the inputs
# --------------------------------------------------------------------------


def dryrun_case(
    directory: Path, scratch: Path, *, np_ranks: int = 4, timeout_s: int = 900
) -> DryRunResult:
    """Run one electronic step in a scratch copy and read back what VASP saw.

    This is the only check that models VASP's own parser.  It is cheap: NELM=1
    stops after the first electronic iteration.

    NSW is capped as well as NELM.  Without that cap a RELAXATION input
    (IBRION=2, NSW=200) runs 200 ionic steps inside the gate -- each step cheap
    because NELM=1, but 200 of them, which turns a two-minute check into an
    hour.  The gate exists to prove VASP can PARSE the input, not to relax it,
    so one ionic step is all it needs.
    """

    if scratch.exists():
        shutil.rmtree(scratch)
    shutil.copytree(directory, scratch)
    incar = scratch / "INCAR"
    text = incar.read_text()
    for tag in ("NELM", "NELMIN", "NSW", "LWAVE", "LCHARG"):
        text = re.sub(rf"^\s*{tag}\s*=.*$", "", text, flags=re.M)
    incar.write_text(
        text + "\nNELM = 1\nNSW = 0\nLWAVE = .FALSE.\nLCHARG = .FALSE.\n"
    )
    for junk in ("WAVECAR", "CHGCAR", "CHG"):
        stale = scratch / junk
        if stale.exists():
            stale.unlink()
    command = (
        "set -e; . /etc/profile.d/modules.sh; ulimit -s unlimited; "
        "module purge >/dev/null 2>&1; module load vasp-cpu >/dev/null 2>&1; "
        f"cd {scratch} && srun -n {np_ranks} vasp_std > dryrun.out 2>&1"
    )
    try:
        subprocess.run(["bash", "-lc", command], timeout=timeout_s, check=False)
    except subprocess.TimeoutExpired:
        return DryRunResult(dir=str(directory), ok=False, problems=["dry run timed out"])
    outcar = scratch / "OUTCAR"
    if not outcar.is_file():
        return DryRunResult(
            dir=str(directory), ok=False, problems=["no OUTCAR produced"],
            scratch=str(scratch),
        )
    text = outcar.read_text(errors="replace")

    def first_int(pattern: str) -> int | None:
        match = re.search(pattern, text)
        return int(match.group(1)) if match else None

    problems = [f"VASP reported: {bad}" for bad in VASP_REJECTIONS if bad in text]
    nions = first_int(r"NIONS\s*=\s*(\d+)")
    if nions is None:
        problems.append(
            "OUTCAR has no NIONS line: VASP did not get as far as reading the cell"
        )
    return DryRunResult(
        dir=str(directory),
        ok=not problems,
        problems=problems,
        scratch=str(scratch),
        nions=nions,
        nkpts=first_int(r"NKPTS\s*=\s*(\d+)"),
        ispin=first_int(r"ISPIN\s*=\s*(\d+)"),
        magmom_echo=first_int(r"MAGMOM.*?NIONS\s*=\s*(\d+)"),
    )


def iter_cases(root: Path) -> list[Path]:
    """Every directory under ``root`` that holds an INCAR.

    Hidden directories are skipped: the dry run writes its scratch copies under
    ``.dryscratch`` inside the same tree, and a second gate run would otherwise
    re-validate (and re-run) its own leftovers.
    """

    return sorted(
        path.parent
        for path in root.rglob("INCAR")
        if not any(part.startswith(".") for part in path.parent.relative_to(root).parts)
    )


def run_gate(
    directories: Sequence[Path], *, scratch: Path | None = None, np_ranks: int = 4
) -> list[GateResult]:
    """Static checks on every directory, then the dry run on those that pass.

    The dry run is skipped for a directory the static layer already rejected:
    there is nothing to learn from asking VASP to parse inputs that are known
    to be malformed, and each run costs a compute allocation.
    """

    results = [validate_dir(d) for d in directories]
    gated: list[GateResult] = []
    for result in results:
        if not result.ok:
            gated.append(result)
            continue
        if scratch is None:
            raise ValueError("run_gate needs a scratch directory to dry-run in")
        case_scratch = scratch / "__".join(Path(result.dir).parts[-2:])
        dry = dryrun_case(Path(result.dir), case_scratch, np_ranks=np_ranks)
        if dry.ok:
            gated.append(result.model_copy(update={"dryrun": dry}))
            continue
        gated.append(
            result.model_copy(
                update={
                    "ok": False,
                    "problems": result.problems
                    + tuple(f"dry run: {problem}" for problem in dry.problems),
                    "dryrun": dry,
                }
            )
        )
    return gated


__all__ = [
    "BOOL_TAGS",
    "DryRunResult",
    "FLOAT_TAGS",
    "GateResult",
    "INCAR_NAMES",
    "INT_TAGS",
    "IncarInfo",
    "KpointsInfo",
    "PoscarInfo",
    "PotcarInfo",
    "VASP_REJECTIONS",
    "dryrun_case",
    "iter_cases",
    "read_incar",
    "read_kpoints",
    "read_poscar",
    "_same_element",
    "parse_repeated_floats",
    "read_potcar",
    "run_gate",
    "validate_dir",
]
