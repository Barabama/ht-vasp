#!/usr/bin/env python
"""Unified VASP input-contract validator.

Why this exists
---------------
Three consecutive incidents in this project had the same shape: a builder wrote
a *template* and the check asked only whether *a key was present*.

  1. POSCAR lattice rows written as literal "{lat:.16f}" strings -> VASP read
     "0 types and 0 ions" and then failed with an unrelated LUSE_VDW message.
  2. INCAR got "ISMEAR = ISMEAR" because a dict was unpacked as a 2-tuple ->
     20 directories "passed" a check that only looked for the substring.
  3. A 16-atom cell got the 2-atom MAGMOM -> VASP refused with
     "Error reading item MAGMOM from file INCAR. Error code was IERR= 5",
     killing 14 cases. The validator had checked POTCAR order, atom counts,
     POSCAR numerics, KPOINTS and ISMEAR -- everything except MAGMOM length.

So the contract here is: **check the VALUE of every tag, against the structure
the directory actually contains**, not the presence of a key.  Every check is
expressed as (a) what was read, (b) what the contract requires, (c) pass/fail.

Two layers
----------
``validate_dir(d)``   static: pure text/structure checks, no VASP needed.
``dryrun_case(d)``    dynamic: copy to scratch, set NELM=1, run one electronic
                      step, and confirm VASP itself echoed NIONS / NKPTS /
                      ISPIN / MAGMOM without an "Error reading item".  This is
                      the only check that can catch a VASP-specific rule the
                      static layer does not model.  It needs a compute node, so
                      it is invoked from a submitted job, not from the CLI here.

CLI
---
    python input_contract.py DIR [DIR ...]              # static checks
    python input_contract.py --root DIR                 # every case under DIR
    python input_contract.py --root DIR --json out.json
Exit status is 1 if any directory fails.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Tags whose value must be an integer / a float.  A tag present with a
# non-literal value (e.g. "ISMEAR = ISMEAR") is a failure, not a warning.
INT_TAGS = ("ISMEAR", "ISPIN", "IBRION", "NSW", "ISIF", "NELM", "NELMIN",
            "NCORE", "KPAR", "NFREE", "ISYM", "ISTART", "ICHARG", "LMAXMIX",
            "NELMDL")
FLOAT_TAGS = ("SIGMA", "ENCUT", "EDIFF", "EDIFFG", "POTIM", "SYMPREC")
BOOL_TAGS = ("LWAVE", "LCHARG", "LREAL", "LOPTICS", "LVTOT", "LDAU", "LASPH")

# A directory may hold several INCAR variants for one case; these are the names
# the project uses for the VASP run itself.
INCAR_NAMES = ("INCAR",)


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()[:16]


# --------------------------------------------------------------------------
# POSCAR
# --------------------------------------------------------------------------
def read_poscar(d: Path) -> dict:
    """Parse just enough POSCAR to know atom counts, species and order."""
    p = d / "POSCAR"
    if not p.is_file():
        return {"error": "no POSCAR"}
    lines = p.read_text().splitlines()
    if len(lines) < 8:
        return {"error": f"POSCAR has only {len(lines)} lines"}
    out: dict = {"lines": len(lines)}
    try:
        out["scale"] = float(lines[1].split()[0])
    except (ValueError, IndexError):
        return {"error": f"POSCAR line 2 is not a scale factor: {lines[1]!r}"}
    lat = []
    for i in (2, 3, 4):
        try:
            v = [float(x) for x in lines[i].split()]
        except ValueError:
            return {"error": f"POSCAR line {i + 1} is not numeric: {lines[i]!r}"}
        if len(v) != 3:
            return {"error": f"POSCAR line {i + 1} has {len(v)} numbers, not 3"}
        lat.append(v)
    out["lattice"] = lat
    # optional "Selective dynamics" line
    idx = 5
    if lines[idx].strip()[:1].lower() == "s":
        idx += 1
    species = lines[idx].split()
    idx += 1
    try:
        counts = [int(x) for x in lines[idx].split()]
    except ValueError:
        return {"error": f"POSCAR count line is not integers: {lines[idx]!r}"}
    if len(counts) != len(species):
        return {"error": f"POSCAR lists {len(species)} species but {len(counts)} counts"}
    out["species"] = species
    out["counts"] = counts
    out["natoms"] = sum(counts)
    out["natoms_declared"] = out["natoms"]
    # coordinate block
    mode = lines[idx + 1].strip()[:1].lower()
    out["coord_mode"] = "Direct" if mode == "d" else ("Cartesian" if mode == "c" else mode)
    body = [ln for ln in lines[idx + 2:] if ln.strip()]
    out["n_coord_lines"] = len(body)
    # per-species MAGMOM expansion must follow POSCAR order
    out["species_order"] = list(dict.fromkeys(species))
    out["per_atom_species"] = [s for s, c in zip(species, counts) for _ in range(c)]
    # volume
    a, b, c = lat
    vol = (a[0] * (b[1] * c[2] - b[2] * c[1])
           - a[1] * (b[0] * c[2] - b[2] * c[0])
           + a[2] * (b[0] * c[1] - b[1] * c[0]))
    out["volume_angstrom3"] = abs(vol) * abs(out["scale"]) ** 3
    out["sha256"] = _sha(p.read_bytes())
    return out


# --------------------------------------------------------------------------
# POTCAR
# --------------------------------------------------------------------------
def read_potcar(d: Path) -> dict:
    p = d / "POTCAR"
    if not p.is_file():
        return {"error": "no POTCAR"}
    txt = p.read_text(errors="replace")
    titles = re.findall(r"TITEL\s*=\s*(\S+)\s+(\w+)", txt)
    zvals = re.findall(r"ZVAL\s*=\s*([\d.]+)", txt)
    enmax = re.findall(r"ENMAX\s*=\s*([\d.]+)", txt)
    return {
        "blocks": len(titles),
        "titles": [t[0] for t in titles],
        "elements": [t[1] for t in titles],
        "zval": [float(z) for z in zvals],
        "enmax": [float(e) for e in enmax],
        "sha256": _sha(p.read_bytes()),
    }


# --------------------------------------------------------------------------
# KPOINTS
# --------------------------------------------------------------------------
def read_kpoints(d: Path) -> dict:
    p = d / "KPOINTS"
    if not p.is_file():
        return {"error": "no KPOINTS"}
    lines = [ln for ln in p.read_text().splitlines() if ln.strip() != ""]
    if len(lines) < 3:
        return {"error": f"KPOINTS has only {len(lines)} non-blank lines"}
    out: dict = {"scheme_line": lines[2].strip(), "sha256": _sha(p.read_bytes())}
    # Automatic mesh: line 4 is the mesh, optional line 5 is the shift.
    if "monkhorst" in lines[2].lower() or "gamma" in lines[2].lower():
        try:
            mesh = [int(x) for x in lines[3].split()]
        except (ValueError, IndexError):
            return {"error": f"KPOINTS mesh line is not integers: {lines[3]!r}"}
        if len(mesh) != 3:
            return {"error": f"KPOINTS mesh has {len(mesh)} entries, not 3"}
        out["mesh"] = mesh
        out["mode"] = "automatic"
        if len(lines) >= 5:
            try:
                out["shift"] = [float(x) for x in lines[4].split()]
            except ValueError:
                return {"error": f"KPOINTS shift line is not numeric: {lines[4]!r}"}
        else:
            out["shift"] = [0.0, 0.0, 0.0]
        if len(out["shift"]) != 3:
            return {"error": f"KPOINTS shift has {len(out['shift'])} entries, not 3"}
        # Gamma inclusion: for a Monkhorst-Pack mesh, Gamma is in the grid iff
        # every n_i is odd (shift 0) or every n_i is even (shift 0.5).  Mixed
        # parity with shift 0.5 is NOT Gamma-centred -- a common silent error.
        n = out["mesh"]
        s = out["shift"]
        all_odd = all(x % 2 == 1 for x in n)
        all_even = all(x % 2 == 0 for x in n)
        if all(x == 0.0 for x in s):
            out["gamma_centered"] = all_odd
            out["gamma_rule"] = "odd mesh + shift 0" if all_odd else "even mesh + shift 0 excludes Gamma"
        elif all(abs(x - 0.5) < 1e-9 for x in s):
            out["gamma_centered"] = all_even
            out["gamma_rule"] = "even mesh + shift 0.5" if all_even else "odd mesh + shift 0.5 excludes Gamma"
        else:
            out["gamma_centered"] = None
            out["gamma_rule"] = f"non-standard shift {s}: Gamma membership not modelled"
    else:
        out["mode"] = "explicit"
        out["mesh"] = None
        out["shift"] = None
        out["gamma_centered"] = None
        out["gamma_rule"] = "explicit list: Gamma membership must be checked in IBZKPT"
    return out


# --------------------------------------------------------------------------
# INCAR
# --------------------------------------------------------------------------
def read_incar(d: Path) -> dict:
    for name in INCAR_NAMES:
        p = d / name
        if p.is_file():
            break
    else:
        return {"error": "no INCAR"}
    raw = p.read_text()
    tags: dict[str, str] = {}
    for ln in raw.splitlines():
        ln = ln.split("#")[0].split("!")[0].strip()
        if not ln or "=" not in ln:
            continue
        k, v = ln.split("=", 1)
        tags[k.strip().upper()] = v.strip()
    return {"tags": tags, "raw": raw, "sha256": _sha(p.read_bytes())}


# --------------------------------------------------------------------------
# the contract
# --------------------------------------------------------------------------
def validate_dir(d: Path, *, expect_natoms: int | None = None,
                 expect_species: list[str] | None = None) -> dict:
    """Return {"ok": bool, "problems": [...], ...} for one VASP directory."""
    problems: list[str] = []
    warns: list[str] = []

    pos = read_poscar(d)
    pot = read_potcar(d)
    kpt = read_kpoints(d)
    inc = read_incar(d)

    for name, obj in (("POSCAR", pos), ("POTCAR", pot),
                      ("KPOINTS", kpt), ("INCAR", inc)):
        if "error" in obj:
            problems.append(f"{name}: {obj['error']}")
    if problems:
        return {"dir": str(d), "ok": False, "problems": problems, "warnings": warns,
                "poscar": pos, "potcar": pot, "kpoints": kpt,
                "incar": {k: v for k, v in inc.items() if k != "raw"}}

    n = pos["natoms"]

    # ---- POSCAR internals -------------------------------------------------
    if pos["n_coord_lines"] != n:
        problems.append(f"POSCAR declares {n} atoms but lists {pos['n_coord_lines']} "
                        f"coordinate lines")
    if expect_natoms is not None and n != expect_natoms:
        problems.append(f"POSCAR has {n} atoms, expected {expect_natoms}")
    if expect_species is not None:
        if [s.upper() for s in pos["species_order"]] != [s.upper() for s in expect_species]:
            problems.append(f"POSCAR species order {pos['species_order']} != expected "
                            f"{expect_species}")
    if any(x == 0 for x in pos["counts"]):
        problems.append(f"POSCAR has a zero species count: {pos['counts']}")
    if abs(pos["scale"]) == 0:
        problems.append("POSCAR scale factor is 0")
    if pos["volume_angstrom3"] <= 0:
        problems.append(f"POSCAR volume is not positive: {pos['volume_angstrom3']}")

    # ---- POTCAR vs POSCAR -------------------------------------------------
    if pot["blocks"] != len(pos["species"]):
        problems.append(f"POTCAR has {pot['blocks']} blocks but POSCAR declares "
                        f"{len(pos['species'])} species")
    if [s.upper() for s in pot["elements"]] != [s.upper() for s in pos["species"]]:
        problems.append(f"POTCAR element order {pot['elements']} != POSCAR order "
                        f"{pos['species']} (VASP matches them positionally)")
    if len(pot["zval"]) != pot["blocks"]:
        problems.append(f"POTCAR has {pot['blocks']} blocks but {len(pot['zval'])} ZVAL values")

    # ---- KPOINTS vs POSCAR ------------------------------------------------
    if kpt["mode"] == "automatic":
        mesh, shift = kpt["mesh"], kpt["shift"]
        if min(mesh) < 1:
            problems.append(f"KPOINTS mesh {mesh} has a non-positive entry")
        if kpt["gamma_centered"] is False and all(abs(x) < 1e-9 for x in shift):
            warns.append(f"KPOINTS mesh {mesh} with shift 0 does NOT include Gamma "
                         f"({kpt['gamma_rule']}); state this explicitly if intended")
        if kpt["gamma_centered"] is None:
            warns.append(f"KPOINTS {kpt['gamma_rule']}")
    else:
        warns.append("KPOINTS is an explicit list; Gamma membership must be read "
                     "from IBZKPT")

    # ---- INCAR value contracts -------------------------------------------
    tags = inc["tags"]
    for t in INT_TAGS:
        if t in tags:
            try:
                int(tags[t])
            except ValueError:
                problems.append(f"INCAR {t} = {tags[t]!r} is not an integer literal")
    for t in FLOAT_TAGS:
        if t in tags:
            try:
                float(tags[t])
            except ValueError:
                problems.append(f"INCAR {t} = {tags[t]!r} is not numeric")
    for t in BOOL_TAGS:
        if t in tags:
            v = tags[t].strip(".").upper()
            if v not in ("TRUE", "FALSE", "T", "F"):
                problems.append(f"INCAR {t} = {tags[t]!r} is not a logical value")

    # MAGMOM: length must be 1 or natoms, and the expansion must match the
    # POSCAR species order.
    def as_int(tag: str, default: int = 0) -> int | None:
        """int value of an INCAR tag, or None if absent/unparseable.

        The INT_TAGS loop above already reports unparseable values; these
        helpers must not raise a second time on the same bad input.
        """
        if tag not in tags:
            return default
        try:
            return int(tags[tag])
        except ValueError:
            return None

    ispin = as_int("ISPIN", 1)
    if "MAGMOM" in tags:
        vals = tags["MAGMOM"].split()
        try:
            fv = [float(v) for v in vals]
        except ValueError:
            problems.append(f"INCAR MAGMOM has a non-numeric value: {tags['MAGMOM']!r}")
            fv = []
        if fv and len(fv) not in (1, n):
            problems.append(f"INCAR MAGMOM has {len(fv)} values but the cell has {n} "
                            f"atoms (VASP accepts only 1 or NIONS)")
        if ispin == 1:
            warns.append("ISPIN=1 with MAGMOM set: MAGMOM is ignored")
    else:
        if ispin == 2:
            warns.append("ISPIN=2 with no MAGMOM: VASP defaults to 1.0 per atom")

    # NCORE / KPAR conflict: both set is a classic VASP parallelisation clash.
    if "NCORE" in tags and "KPAR" in tags:
        problems.append("INCAR sets both NCORE and KPAR; VASP documents them as "
                        "mutually exclusive parallelisation controls")

    # ICHARG / LCHARG consistency with the restart story
    icharg = as_int("ICHARG", 0)
    if icharg in (1, 11) and "CHGCAR" not in {f.name for f in d.iterdir()}:
        problems.append(f"ICHARG={tags['ICHARG']} needs a CHGCAR in the directory; none present")
    if icharg == 1 and as_int("ISTART", 0) != 1:
        warns.append(f"ICHARG=1 with ISTART={tags.get('ISTART', '0')}; VASP resets ISTART to 1")

    # IBRION sanity: 8 = DFPT, -1 = static, 6 = finite difference
    ib = as_int("IBRION", -1)
    if ib is not None and ib not in (-1, 0, 1, 2, 3, 5, 6, 7, 8):
        warns.append(f"IBRION={ib} is unusual")

    return {
        "dir": str(d),
        "ok": not problems,
        "problems": problems,
        "warnings": warns,
        "natoms": n,
        "species_order": pos["species_order"],
        "volume_angstrom3": pos["volume_angstrom3"],
        "poscar_sha256": pos["sha256"],
        "potcar_sha256": pot["sha256"],
        "incar_sha256": inc["sha256"],
        "kpoints_sha256": kpt["sha256"],
        "potcar_zval": pot["zval"],
        "mesh": kpt["mesh"],
        "shift": kpt["shift"],
        "gamma_centered": kpt["gamma_centered"],
        "ismear": tags.get("ISMEAR"),
        "sigma": tags.get("SIGMA"),
    }


# --------------------------------------------------------------------------
# dynamic gate: let VASP itself parse the inputs
# --------------------------------------------------------------------------
def dryrun_case(d: Path, scratch: Path, *, np_ranks: int = 4,
                timeout_s: int = 900) -> dict:
    """Run one electronic step in a scratch copy and read back what VASP saw.

    This is the only check that models VASP's own parser, which is where the
    MAGMOM-length and POSCAR-layout failures came from.  It is cheap: NELM=1
    stops after the first electronic iteration.

    NSW is capped as well as NELM.  Without that cap a RELAXATION input
    (IBRION=2, NSW=200) runs 200 ionic steps inside the gate -- each step cheap
    because NELM=1, but 200 of them, which turns a two-minute check into an
    hour.  The gate exists to prove VASP can PARSE the input, not to relax it,
    so one ionic step is all it needs.
    """
    if scratch.exists():
        shutil.rmtree(scratch)
    shutil.copytree(d, scratch)
    # force a single iteration and suppress large outputs
    incar = scratch / "INCAR"
    txt = incar.read_text()
    txt = re.sub(r"^\s*NELM\s*=.*$", "", txt, flags=re.M)
    txt = re.sub(r"^\s*NELMIN\s*=.*$", "", txt, flags=re.M)
    # cap the IONIC loop too -- see the docstring
    txt = re.sub(r"^\s*NSW\s*=.*$", "", txt, flags=re.M)
    txt = re.sub(r"^\s*LWAVE\s*=.*$", "", txt, flags=re.M)
    txt = re.sub(r"^\s*LCHARG\s*=.*$", "", txt, flags=re.M)
    incar.write_text(txt + "\nNELM = 1\nNSW = 0\n"
                            "LWAVE = .FALSE.\nLCHARG = .FALSE.\n")
    for junk in ("WAVECAR", "CHGCAR", "CHG"):
        f = scratch / junk
        if f.exists():
            f.unlink()
    cmd = ("set -e; . /etc/profile.d/modules.sh; ulimit -s unlimited; "
           "module purge >/dev/null 2>&1; module load vasp-cpu >/dev/null 2>&1; "
           f"cd {scratch} && srun -n {np_ranks} vasp_std > dryrun.out 2>&1")
    try:
        subprocess.run(["bash", "-lc", cmd], timeout=timeout_s, check=False)
    except subprocess.TimeoutExpired:
        return {"dir": str(d), "ok": False, "problems": ["dry run timed out"]}
    out = scratch / "OUTCAR"
    res: dict = {"dir": str(d), "ok": False, "problems": [], "scratch": str(scratch)}
    if not out.is_file():
        res["problems"].append("no OUTCAR produced")
        return res
    txt = out.read_text(errors="replace")
    res["nions"] = (int(m.group(1)) if (m := re.search(r"NIONS\s*=\s*(\d+)", txt)) else None)
    res["nkpts"] = (int(m.group(1)) if (m := re.search(r"NKPTS\s*=\s*(\d+)", txt)) else None)
    res["ispin"] = (int(m.group(1)) if (m := re.search(r"ISPIN\s*=\s*(\d+)", txt)) else None)
    res["magmom_echo"] = (int(m.group(1)) if (m := re.search(r"MAGMOM.*?NIONS\s*=\s*(\d+)", txt)) else None)
    for bad in ("Error reading item", "I REFUSE TO CONTINUE",
                "POSCAR found : 0 types", "not enough POTCAR"):
        if bad in txt:
            res["problems"].append(f"VASP reported: {bad}")
    if res["nions"] is None:
        res["problems"].append("OUTCAR has no NIONS line: VASP did not get as far as reading the cell")
    res["ok"] = not res["problems"]
    return res


def iter_cases(root: Path) -> list[Path]:
    """Every directory under root that holds an INCAR.

    Hidden directories are skipped: ``--dryrun`` writes its scratch copies under
    ``.dryscratch`` inside the same tree, and a second gate run would otherwise
    re-validate (and re-run) its own leftovers.
    """
    return sorted(p.parent for p in root.rglob("INCAR")
                  if not any(part.startswith(".") for part in p.parent.relative_to(root).parts))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="*", type=Path)
    ap.add_argument("--root", type=Path, default=None,
                    help="validate every directory under here that has an INCAR")
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--dryrun", action="store_true",
                    help="also let VASP parse each case (needs a compute allocation; "
                         "each case is copied to --scratch and run with NELM=1)")
    ap.add_argument("--scratch", type=Path, default=None,
                    help="scratch root for --dryrun copies")
    ap.add_argument("--ranks", type=int, default=4)
    args = ap.parse_args()

    targets: list[Path] = list(args.dirs)
    if args.root:
        targets += iter_cases(args.root)
    if not targets:
        ap.error("give directories or --root")

    results = [validate_dir(d) for d in targets]
    if args.dryrun:
        if args.scratch is None:
            ap.error("--dryrun needs --scratch")
        args.scratch.mkdir(parents=True, exist_ok=True)
        for r in results:
            if not r["ok"]:
                r["dryrun"] = {"skipped": "static checks failed"}
                continue
            d = Path(r["dir"])
            # unique scratch name: basenames repeat across endpoints
            # (d16_k16_s1 exists for several representatives)
            key = "__".join(d.parts[-2:]) if len(d.parts) >= 2 else d.name
            dr = dryrun_case(d, args.scratch / key, np_ranks=args.ranks)
            r["dryrun"] = dr
            if not dr["ok"]:
                r["ok"] = False
                r["problems"] = list(r["problems"]) + [
                    f"dry run: {p}" for p in dr["problems"]]

    nbad = sum(1 for r in results if not r["ok"])
    if not args.quiet:
        for r in results:
            mark = "OK  " if r["ok"] else "FAIL"
            extra = ""
            if r.get("mesh"):
                extra = (f" mesh={r['mesh']} shift={r['shift']} "
                         f"gamma={r['gamma_centered']} ISMEAR={r.get('ismear')}")
            print(f"{mark} {r['dir']}{extra}")
            for p in r["problems"]:
                print(f"       !! {p}")
            for w in r.get("warnings", []):
                print(f"       ~  {w}")
            if "dryrun" in r and isinstance(r["dryrun"], dict) and "nions" in r["dryrun"]:
                dr = r["dryrun"]
                print(f"       VASP read: NIONS={dr['nions']} NKPTS={dr['nkpts']} "
                      f"ISPIN={dr['ispin']}")
        print(f"\n{len(results) - nbad}/{len(results)} directories pass the contract")
    if args.json:
        args.json.write_text(json.dumps(results, indent=2) + "\n")
        print(f"wrote {args.json}")
    if nbad:
        sys.exit(1)


if __name__ == "__main__":
    main()
