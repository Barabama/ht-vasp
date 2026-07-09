#!/usr/bin/env python3
"""Calculate primitive-cell lattice constants with ISIF=8 for each endmember.

Usage:
    python calc_isif8.py                         # list all structures
    python calc_isif8.py --run SER-Fe            # run locally
    python calc_isif8.py --slurm                 # submit all to Slurm
    python calc_isif8.py --slurm --dry           # dry-run (print sbatch commands)
    python calc_isif8.py --gen-only SER-Co       # only write input files
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import time
from pathlib import Path
from datetime import datetime

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s"
)
log = logging.getLogger(__name__)

WORK_DIR = Path(__file__).parent.resolve()
POSCAR_DIR = WORK_DIR / "data" / "poscars"
OUTPUT_DIR = WORK_DIR / "data" / "isif8"

# Structures to process (all available POSCARs)
# Auto-detected from data/poscars/*.vasp

# VASP command (CPU partition)
VASP_CMD = (
    "/bin/bash -c '"
    ". /etc/profile.d/modules.sh && "
    "module load vasp-cpu && "
    "srun vasp_std"
    "'"
)

ISIF8_INCAR = {
    "SYSTEM": "ISIF=8 lattice constant scan",
    "ENCUT": 500,
    "ISTART": 0,
    "ICHARG": 2,
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "ALGO": "Normal",
    "NELM": 120,
    "NELMIN": 6,
    "NELMDL": -6,
    "IBRION": 2,
    "ISIF": 8,
    "NSW": 30,
    "POTIM": 0.2,
    "EDIFF": 1e-6,
    "EDIFFG": -0.02,
    "ISYM": 0,
    "LREAL": "Auto",
    "PREC": "Normal",
    "SYMPREC": 1e-5,
    "LWAVE": False,
    "LCHARG": False,
    "LORBIT": 10,
    "LAECHG": False,
    "GGA": "PE",
    "KPAR": 2,
    "NCORE": 4,
    "AMIX": 0.1,
    "BMIX": 0.0001,
    "AMIX_MAG": 0.4,
    "BMIX_MAG": 0.0001,
}

# Structures to process
# Combines data/poscars/ and qha_AlNbTiV/data/poscars/
POSCAR_DIRS = [
    POSCAR_DIR,
    Path("qha_AlNbTiV/data/poscars"),
]
ALL_NAMES = set()
for d in POSCAR_DIRS:
    if d.exists():
        for f in d.glob("*.vasp"):
            ALL_NAMES.add(f.stem)

PREFERRED = {
    # Al-Nb-Ti-V BCC (current project)
    "SER-Al", "SER-Nb", "SER-Ti", "SER-V",
    "BCC-Al-Al", "BCC-Al-Nb", "BCC-Al-Ti", "BCC-Al-V",
    "BCC-Nb-Nb", "BCC-Nb-Ti", "BCC-Nb-V",
    "BCC-Ti-Ti", "BCC-Ti-V",
    "BCC-V-V",
    # Co-Fe BCC/FCC (OJ extension)
    "SER-Co", "SER-Fe",
    "BCC-Co-Co", "BCC-Fe-Fe", "BCC-Co-Fe",
    "FCC-Co-Co", "FCC-Fe-Fe", "FCC-Co-Fe", "FCC-Fe-Ni",
}
STRUCT_NAMES = sorted(ALL_NAMES & PREFERRED)


def read_poscar(name: str) -> str:
    """Read POSCAR file content from any known source directory."""
    for d in POSCAR_DIRS:
        path = d / f"{name}.vasp"
        if path.exists():
            return path.read_text()
    raise FileNotFoundError(f"POSCAR not found for {name} in {POSCAR_DIRS}")


def write_vasp_inputs(name: str, poscar_text: str, incar: dict, out_dir: Path):
    """Write INCAR, POSCAR, KPOINTS to out_dir.

    POTCAR is NOT written here — it's handled per-element setup.
    Returns: Path to the output directory.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # INCAR
    incar_lines = []
    for k, v in incar.items():
        if isinstance(v, bool):
            incar_lines.append(f"{k} = {'.TRUE.' if v else '.FALSE.'}")
        elif isinstance(v, str):
            incar_lines.append(f"{k} = {v}")
        else:
            incar_lines.append(f"{k} = {v}")
    (out_dir / "INCAR").write_text("\n".join(incar_lines) + "\n")

    # POSCAR
    (out_dir / "POSCAR").write_text(poscar_text)

    # KPOINTS (Gamma-centered, 4x4x4 minimum, scale by cell size)
    (out_dir / "KPOINTS").write_text("Auto\n0\nGamma\n4 4 4\n0 0 0\n")

    # PBS script
    pbs_lines = [
        "#!/bin/bash",
        "#SBATCH --job-name=isif8_" + name,
        "#SBATCH --output=isif8_" + name + ".log",
        "#SBATCH --ntasks=32",
        "#SBATCH --mem=20G",
        "#SBATCH --partition=partCPU",
        "#SBATCH --time=4:00:00",
        "",
        ". /etc/profile.d/modules.sh",
        "ulimit -s unlimited",
        "module purge",
        "module load vasp-cpu",
        "",
        "srun vasp_std",
    ]
    (out_dir / "run_vasp.slurm").write_text("\n".join(pbs_lines) + "\n")

    return out_dir


def extract_lattice(out_dir: Path) -> dict | None:
    """Extract final lattice constants from CONTCAR or OUTCAR."""
    # Try CONTCAR first (relaxed structure)
    contcar = out_dir / "CONTCAR"
    outcar = out_dir / "OUTCAR"

    if contcar.exists():
        from pymatgen.core import Structure
        try:
            struct = Structure.from_file(contcar)
            abc = struct.lattice.abc
            return {
                "a": abc[0],
                "b": abc[1],
                "c": abc[2],
                "alpha": struct.lattice.alpha,
                "beta": struct.lattice.beta,
                "gamma": struct.lattice.gamma,
                "volume": struct.lattice.volume,
                "natoms": len(struct),
                "source": "CONTCAR",
            }
        except Exception as e:
            log.warning(f"  CONTCAR parse failed: {e}")

    # Fallback: OUTCAR last relaxation step
    if outcar.exists():
        import re
        text = outcar.read_text()
        # Look for last lattice constants
        matches = re.findall(
            r"length of vectors\s*\n\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\n"
            r"\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\n"
            r"\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)",
            text,
        )
        if matches:
            last = matches[-1]
            return {
                "a": float(last[0]),
                "b": float(last[4]),
                "c": float(last[8]),
                "volume": None,
                "source": "OUTCAR(last)",
            }

        # Or total energy
        eng_match = re.findall(
            r"energy without entropy\s*=\s*([-\d.]+)", text
        )
        if eng_match:
            return {"energy_eV": float(eng_match[-1]), "source": "OUTCAR(energy only)"}

    return None


def run_single(name: str, out_dir: Path, force: bool = False) -> dict:
    """Run ISIF=8 calculation locally."""
    json_path = out_dir / "result.json"

    if json_path.exists() and not force:
        data = json.loads(json_path.read_text())
        if data.get("a"):
            log.info(f"  ⏭ {name} already done")
            return data

    out_dir.mkdir(parents=True, exist_ok=True)

    poscar_text = read_poscar(name)
    write_vasp_inputs(name, poscar_text, ISIF8_INCAR, out_dir)

    # Check if we need POTCAR
    potcar_path = out_dir / "POTCAR"
    if not potcar_path.exists():
        # Generate POTCAR from POSCAR species
        from pymatgen.core import Structure
        from pymatgen.io.vasp import Potcar
        struct = Structure.from_str(poscar_text, fmt="poscar")
        symbols = [str(s.specie) for s in struct]
        unique_symbols = list(set(symbols))
        try:
            potcar = Potcar(symbols, functional="PBE_64")
            potcar.write_file(str(potcar_path))
        except Exception as e:
            log.error(f"  Cannot generate POTCAR: {e}")
            return {"error": str(e)}

    # Copy POTCAR from a reference directory if generation failed
    if not potcar_path.exists():
        ref_potcar = Path("qha_AlNbTiV/data/SER-Al/10-phonon_static_eos_deformation_1/POTCAR.gz")
        if ref_potcar.exists():
            import gzip
            with gzip.open(ref_potcar, "rb") as f_in, open(potcar_path, "wb") as f_out:
                f_out.write(f_in.read())
            log.warning("  Copied reference POTCAR (may not match species!)")

    log.info(f"  ▶ Running VASP for {name} in {out_dir}")
    os.chdir(out_dir)
    result = subprocess.run(
        ["srun", "vasp_std"],
        capture_output=True, text=True, timeout=3600 * 4,
    )
    os.chdir(WORK_DIR)

    if result.returncode != 0:
        log.error(f"  ❌ VASP failed for {name}")
        return {"error": result.stderr[:500]}

    data = extract_lattice(out_dir) or {}
    json_path.write_text(json.dumps(data, indent=2, default=str))
    return data


def submit_slurm(name: str, out_dir: Path, dry_run: bool = False):
    """Write inputs and submit to Slurm."""
    poscar_text = read_poscar(name)
    write_vasp_inputs(name, poscar_text, ISIF8_INCAR, out_dir)

    # Generate POTCAR
    from pymatgen.core import Structure
    from pymatgen.io.vasp import Potcar
    struct = Structure.from_str(poscar_text, fmt="poscar")
    symbols = [str(s.specie) for s in struct]
    unique_symbols = list(set(symbols))

    potcar_path = out_dir / "POTCAR"
    if not potcar_path.exists():
        try:
            potcar = Potcar(unique_symbols, functional="PBE_64")
            potcar.write_file(str(potcar_path))
        except Exception as e:
            log.warning(f"  POTCAR generation failed: {e}")

    sbatch_path = out_dir / "run_vasp.slurm"
    log.info(f"  {'[DRY]' if dry_run else ''} Submit: sbatch {sbatch_path}")

    if not dry_run:
        import subprocess as sp
        result = sp.run(
            ["sbatch", str(sbatch_path)],
            capture_output=True, text=True, cwd=out_dir,
        )
        if result.returncode == 0:
            for line in result.stdout.splitlines():
                if "Submitted batch job" in line:
                    job_id = line.split()[-1]
                    log.info(f"  ✅ {name}: job {job_id}")
                    return job_id
        log.error(f"  ❌ submission failed: {result.stderr}")
        return None
    return None


def check_results():
    """Scan output dir for completed calculations and print table."""
    results = {}
    for d in sorted(OUTPUT_DIR.iterdir()):
        if not d.is_dir():
            continue
        json_path = d / "result.json"
        name = d.name
        if json_path.exists():
            data = json.loads(json_path.read_text())
            results[name] = data
        else:
            # Check if VASP finished
            outcar = d / "OUTCAR"
            if outcar.exists():
                data = extract_lattice(d)
                if data:
                    json_path.write_text(json.dumps(data, indent=2, default=str))
                    results[name] = data

    if not results:
        print("  No completed results found.")
        print(f"  Output dir: {OUTPUT_DIR}")
        return

    print(f"\n{'Name':25s} {'a (A)':>10} {'b (A)':>10} {'c (A)':>10} {'V (A^3)':>12} {'n_atoms':>8} {'source':>15}")
    print("-" * 92)
    for name in sorted(results):
        data = results[name]
        a = data.get("a", "?")
        b = data.get("b", "?")
        c = data.get("c", "?")
        vol = f"{data['volume']:.3f}" if data.get("volume") else "?"
        nat = data.get("natoms", "?")
        src = data.get("source", "?")
        energy = f"  E={data['energy_eV']:.4f}eV" if data.get("energy_eV") else ""
        print(f"{name:25s} {str(a):>10} {str(b):>10} {str(c):>10} {str(vol):>12} {str(nat):>8} {str(src):>15}{energy}")


def main():
    p = argparse.ArgumentParser(description="ISIF=8 lattice constant calculator")
    p.add_argument("--run", type=str, default=None, help="Run single structure locally")
    p.add_argument("--slurm", action="store_true", help="Submit all to Slurm")
    p.add_argument("--dry", action="store_true", help="Dry-run (with --slurm)")
    p.add_argument("--force", action="store_true", help="Re-run even if done")
    p.add_argument("--results", action="store_true", help="Check completed results")
    p.add_argument("--gen-only", type=str, default=None, help="Only write input files")
    p.add_argument("--name", type=str, default=None, help="Single structure name (with --slurm)")
    args = p.parse_args()

    os.chdir(WORK_DIR)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.results:
        check_results()
        return

    if args.gen_only:
        name = args.gen_only
        out_dir = OUTPUT_DIR / name
        poscar_text = read_poscar(name)
        write_vasp_inputs(name, poscar_text, ISIF8_INCAR, out_dir)
        log.info(f"  Written inputs to {out_dir}")
        return

    if args.run:
        name = args.run
        out_dir = OUTPUT_DIR / name
        run_single(name, out_dir, force=args.force)
        return

    if args.slurm:
        names = [args.name] if args.name else STRUCT_NAMES
        for name in names:
            out_dir = OUTPUT_DIR / name
            result_json = out_dir / "result.json"

            if result_json.exists() and not args.force:
                data = json.loads(result_json.read_text())
                if data.get("a"):
                    log.info(f"  ⏭ {name} already done")
                    continue

            submit_slurm(name, out_dir, dry_run=args.dry)
            time.sleep(0.3)  # avoid hammering sbatch
        return

    # Default: list structures
    print(f"Available structures ({len(STRUCT_NAMES)}):")
    for n in STRUCT_NAMES:
        out_dir = OUTPUT_DIR / n
        result_json = out_dir / "result.json"
        status = "✅" if result_json.exists() and json.loads(result_json.read_text()).get("a") else "⏳"
        print(f"  {status} {n}")

    print(f"\nCommands:")
    print(f"  python {__file__} --gen-only SER-Fe    # write inputs only")
    print(f"  python {__file__} --run SER-Fe         # run locally")
    print(f"  python {__file__} --slurm --dry        # dry-run")
    print(f"  python {__file__} --slurm              # submit all to Slurm")
    print(f"  python {__file__} --results            # check results")


if __name__ == "__main__":
    main()
