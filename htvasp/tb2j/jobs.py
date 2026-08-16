"""TB2J Jobs - Job-decorated function for the TB2J exchange solve step."""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from jobflow import job
from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Vasprun

log = logging.getLogger(__name__)

# Number pattern for the exchange.out numeric fields. Handles integers,
# decimals, trailing-dot notation ("0.") and scientific notation ("1.2e-05").
# TB2J could emit any of these after a version change, and the pair parser must
# keep matching the whole line instead of silently dropping the pair.
_NUM_PATTERN = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
# Matches numbers in tensor rows, including trailing-dot notation like "0.   ".
_NUM_RE = re.compile(_NUM_PATTERN)

# Line like: "   Co1   Ni1   ( -1,  -1,  -1)  4.8244   (-1.406, -1.406, -1.406)  2.435 "
# The jiso/vector/distance numeric groups share _NUM_PATTERN with _NUM_RE so
# trailing-dot / scientific-notation values match the whole line. The cell
# image R is always an integer (-?\d+).
_PAIR_RE = re.compile(
    r"^\s*(?P<i>[A-Z][a-z]?\d+)\s+(?P<j>[A-Z][a-z]?\d+)\s+"
    r"\(\s*(?P<r0>-?\d+)\s*,\s*(?P<r1>-?\d+)\s*,\s*(?P<r2>-?\d+)\s*\)\s+"
    rf"(?P<jiso>{_NUM_PATTERN})\s+"
    rf"\(\s*(?P<v0>{_NUM_PATTERN})\s*,\s*(?P<v1>{_NUM_PATTERN})\s*,\s*(?P<v2>{_NUM_PATTERN})\s*\)\s+"
    rf"(?P<dist>{_NUM_PATTERN})\s*$"
)


def _strip_host(prev_dir: str | Path) -> Path:
    """Strip a jobflow ``host:path`` prefix from a directory reference.

    jobflow's local runner stores job ``dir_name`` values host-qualified (e.g.
    ``429pro:/nfs_ssd/tmp/SER-Co-cf44920f/3-w90_static``). Downstream jobs that
    receive ``job.output.dir_name`` must drop the ``host:`` prefix before
    touching the filesystem. Plain paths are returned unchanged.

    Args:
        prev_dir: Possibly host-qualified path.

    Returns:
        Path without any ``host:`` prefix.
    """
    s = str(prev_dir)
    head = s.split("/", 1)[0]
    if ":" in head:
        s = s.split(":", 1)[1]
    return Path(s)


def _ensure_file(prev_dir: Path, name: str) -> Path:
    """Return ``prev_dir/name``, decompressing ``name.gz`` if necessary.

    atomate2's ``BaseVaspMaker`` gzips the VASP input/output files of every
    completed job (``SETTINGS.VASP_ZIP_FILES`` = "atomate", covering e.g.
    vasprun.xml, POSCAR, OUTCAR) *before* the downstream TB2J solve job runs.
    Restore the plain file so ``Vasprun`` / TB2J's ``wann2J`` can read it.

    Args:
        prev_dir: The (host-prefix-stripped) Wannier90 output directory.
        name: File name to ensure exists uncompressed.

    Returns:
        Path to the plain file.
    """
    plain = prev_dir / name
    if plain.exists():
        return plain
    gz = prev_dir / f"{name}.gz"
    if gz.exists():
        import gzip

        with gzip.open(gz, "rb") as fin, open(plain, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        log.debug(f"Decompressed {gz.name} -> {name} for TB2J solve")
    return plain


def get_efermi(prev_dir: str | Path) -> float | None:
    """Extract the Fermi energy from a VASP ``vasprun.xml``.

    Args:
        prev_dir: Directory containing the SCF + Wannier90 vasprun.xml.

    Returns:
        Fermi energy in eV, or None if it could not be read.
    """
    prev_dir = _strip_host(prev_dir)
    vasprun_path = _ensure_file(prev_dir, "vasprun.xml")
    if not vasprun_path.exists():
        log.warning(f"vasprun.xml not found in {prev_dir}")
        return None
    try:
        vasprun = Vasprun(vasprun_path, parse_dos=True, parse_potcar_file=False)
        return vasprun.efermi
    except Exception as exc:  # noqa: BLE001 - report and fall back
        log.warning(f"Failed to read efermi from {vasprun_path}: {exc}")
        return None


def run_wann2j(
    prev_dir: str | Path,
    efermi: float,
    elements: list[str],
    kmesh: tuple[int, int, int] = (9, 9, 9),
    output_path: str | Path | None = None,
) -> subprocess.CompletedProcess:
    """Run TB2J's ``python -m TB2J.scripts.wann2J``.

    Args:
        prev_dir: Directory with the Wannier90 output (wannier90.1_hr.dat, ...).
        efermi: Fermi energy in eV.
        elements: Magnetic elements (space separated for --elements).
        kmesh: Monkhorst-Pack / Gamma grid for the real-space interpolation.
        output_path: Directory where exchange.out is written (default: cwd).

    Returns:
        CompletedProcess from subprocess.run.
    """
    prev_dir = _strip_host(prev_dir).resolve()
    output_path = Path(output_path) if output_path is not None else Path.cwd()
    output_path.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        "-m",
        "TB2J.scripts.wann2J",
        "--path",
        str(prev_dir),
        "--posfile",
        str(prev_dir / "POSCAR"),
        "--efermi",
        f"{efermi:.8f}",
        "--kmesh",
        *(str(k) for k in kmesh),
        "--elements",
        *elements,
        "--prefix_up",
        "wannier90.1",
        "--prefix_down",
        "wannier90.2",
        "--output_path",
        str(output_path),
    ]
    log.info(f"Running TB2J: {' '.join(cmd)}")
    return subprocess.run(cmd, capture_output=True, text=True)


def parse_exchange_out(path: str | Path) -> dict[str, Any]:
    """Parse a TB2J ``exchange.out`` file into structured J data.

    Args:
        path: Path to exchange.out.

    Returns:
        Dict with ``j_pairs`` (list of {i, j, R, J_iso_meV, vector, distance}),
        ``j_tensors`` (list of 3x3 meV tensors, parallel to j_pairs),
        ``atom_moments`` (per-atom wannier charge/magmom), ``cell`` and
        ``warnings`` (non-fatal messages, e.g. J_iso markers with no matching
        pair line).
    """
    text = Path(path).read_text(encoding="utf-8")
    lines = text.splitlines()

    j_pairs: list[dict[str, Any]] = []
    j_tensors: list[list[list[float]]] = []
    atom_moments: list[dict[str, Any]] = []
    cell: list[list[float]] | None = None

    # Line index of the last successfully parsed pair line. A "J_iso:" marker
    # whose immediately preceding line was NOT such a pair means the pair line
    # changed format and was silently skipped - count those so the caller can
    # warn instead of losing J pairs unnoticed.
    last_pair_line: int | None = None
    dropped_jiso: int = 0

    i, n = 0, len(lines)
    while i < n:
        line = lines[i]

        if line.startswith("Cell (Angstrom):"):
            cell = [list(map(float, lines[i + 1 + k].split())) for k in range(3)]
            i += 4
            continue

        if line.startswith("Atom number") and "w_charge" in line:
            i += 1
            while i < n:
                tok = lines[i].strip()
                if not tok or tok.startswith("Total"):
                    break
                parts = tok.split()
                if len(parts) == 6:
                    atom_moments.append(
                        {
                            "atom": parts[0],
                            "xyz": [float(x) for x in parts[1:4]],
                            "w_charge": float(parts[4]),
                            "w_magmom": float(parts[5]),
                        }
                    )
                i += 1
            continue

        if line.startswith("J_iso:") and last_pair_line != i - 1:
            dropped_jiso += 1

        m = _PAIR_RE.match(line)
        if m:
            last_pair_line = i
            d = m.groupdict()
            tensor: list[list[float]] = []
            k = i + 1
            while k < n and len(tensor) < 3:
                if "[" in lines[k]:
                    nums = _NUM_RE.findall(lines[k])
                    if nums:
                        tensor.append([float(x) for x in nums])
                k += 1
            j_pairs.append(
                {
                    "i": d["i"],
                    "j": d["j"],
                    "R": [int(d["r0"]), int(d["r1"]), int(d["r2"])],
                    "J_iso_meV": float(d["jiso"]),
                    "vector": [float(d["v0"]), float(d["v1"]), float(d["v2"])],
                    "distance": float(d["dist"]),
                }
            )
            j_tensors.append(tensor)
        i += 1

    warnings: list[str] = []
    if dropped_jiso:
        warnings.append(
            f"{dropped_jiso} J_iso marker(s) had no matching pair line - "
            "J pair(s) were silently dropped (e.g. trailing-dot or "
            "scientific-notation numbers in exchange.out)"
        )

    return {
        "exchange_out_path": str(Path(path).resolve()),
        "cell": cell,
        "atom_moments": atom_moments,
        "j_pairs": j_pairs,
        "j_tensors": j_tensors,
        "num_j_pairs": len(j_pairs),
        "warnings": warnings,
    }


@job
def tb2j_solve(
    prev_dir: str | Path,
    elements: list[str],
    kmesh: tuple[int, int, int] = (9, 9, 9),
    num_wann: int | None = None,
) -> dict[str, Any]:
    """Solve for exchange parameters J using TB2J (wann2J).

    Reads the Fermi energy from ``prev_dir``/vasprun.xml, runs
    ``python -m TB2J.scripts.wann2J`` with the collinear Wannier90 output, and
    parses the resulting ``exchange.out`` into structured J data.

    Args:
        prev_dir: Directory with the SCF + Wannier90 VASP output.
        elements: Magnetic elements (passed to TB2J --elements).
        kmesh: kmesh for the TB2J real-space interpolation.
        num_wann: Number of Wannier orbitals (for provenance; optional).

    Returns:
        Dict with exchange_out_path, efermi, num_atoms, elements, kmesh,
        num_wann, j_pairs, j_tensors, atom_moments and warnings.
    """
    prev_dir = _strip_host(prev_dir)
    warnings: list[str] = []

    efermi = get_efermi(prev_dir)
    if efermi is None:
        raise RuntimeError(f"Could not extract Fermi energy from {prev_dir}")

    num_atoms = None
    poscar_path = _ensure_file(prev_dir, "POSCAR")
    if poscar_path.exists():
        num_atoms = len(Structure.from_file(poscar_path))

    result = run_wann2j(prev_dir, efermi, elements, kmesh=kmesh, output_path=Path.cwd())
    if result.returncode != 0:
        raise RuntimeError(
            f"TB2J wann2J failed (returncode={result.returncode}): "
            f"{result.stderr.strip()}\nstdout:\n{result.stdout.strip()}"
        )

    exchange_out_path = Path.cwd() / "exchange.out"
    if not exchange_out_path.exists():
        raise RuntimeError(f"TB2J finished but no exchange.out in {Path.cwd()}")

    parsed = parse_exchange_out(exchange_out_path)
    warnings.extend(parsed.get("warnings", []))
    if parsed["num_j_pairs"] == 0:
        warnings.append("No J pairs parsed from exchange.out")

    if num_atoms is None and parsed["atom_moments"]:
        num_atoms = len(parsed["atom_moments"])

    return {
        "exchange_out_path": parsed["exchange_out_path"],
        "efermi": efermi,
        "num_atoms": num_atoms,
        "elements": list(elements),
        "kmesh": list(kmesh),
        "num_wann": num_wann,
        "j_pairs": parsed["j_pairs"],
        "j_tensors": parsed["j_tensors"],
        "atom_moments": parsed["atom_moments"],
        "num_j_pairs": parsed["num_j_pairs"],
        "warnings": warnings,
    }
