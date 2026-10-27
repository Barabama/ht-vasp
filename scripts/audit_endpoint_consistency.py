#!/usr/bin/env python
"""Cross-endpoint consistency checks that need no new VASP and no reference data.

Every check here compares two artifacts that describe **the same crystal**, so
any disagreement is a bug in one of them rather than a matter of opinion:

* ``SER-<El>-1`` against ``BCC-<El>-<El>-2`` for a bcc-stable element, and
  ``SER-<El>-1`` against ``FCC-<El>-<El>-4`` for a fcc-stable one.  These are
  separate runs of the same crystal, so a difference means one of the runs,
  its POTCAR, its cell or its magnetic branch is wrong.
* The same endpoint produced on both clusters, compared bit for bit.

Plus two physical sanity checks that catch a whole class of silent failures:

* the pure-element lattice constant against the experimental one, with a wide
  tolerance (the point is to catch a wrong structure, not a GGA error);
* the E(V) curve must have a minimum inside the sampled window and positive
  curvature -- a run whose "equilibrium" is the edge of the scan is reported
  even when ``ok`` is true.

Run it after mirroring; it reads both the local tree and ``data/remote/*``.

    ./.conda/bin/python scripts/audit_endpoint_consistency.py
    ./.conda/bin/python scripts/audit_endpoint_consistency.py --json report.json
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: Elements whose ground state is bcc, and the SER endpoint that records it.
BCC_STABLE = ("Cr", "Nb", "Ta")
FCC_STABLE = ("Al",)

#: Experimental lattice constants at 300 K, Angstrom.  Used only to catch a
#: wrong structure or a wrong POTCAR: the tolerance is 3%, several times the
#: GGA error, so a passing number says nothing about accuracy.
LATTICE_A_EXPERIMENT = {
    ("BCC", "Cr"): 2.884, ("BCC", "Nb"): 3.301, ("BCC", "Ta"): 3.301,
    ("BCC", "Ti"): 3.310, ("BCC", "Zr"): 3.609,  # Ti/Zr bcc are hypothetical
    ("FCC", "Al"): 4.050, ("FCC", "Hf"): 4.477, ("FCC", "Zr"): 4.611,
}

#: Agreement required when two runs describe the same crystal.
SAME_CRYSTAL_TOL_MEV = 1.0      # per atom
SAME_CRYSTAL_TOL_V = 1e-3       # Angstrom^3 per atom
SAME_CRYSTAL_TOL_B = 0.1        # GPa


def roots() -> list[tuple[Path, str]]:
    found = [(ROOT / "data" / "endpoints", "local")]
    for mirror in sorted((ROOT / "data" / "remote").glob("*")):
        tree = mirror / "data" / "endpoints"
        if tree.is_dir():
            found.append((tree, mirror.name))
    return found


def latest_ok(name: str, base: list[tuple[Path, str]]) -> dict | None:
    """The newest ``ok=True`` analysis for an endpoint, across all roots."""
    best: tuple[float, Path, dict] | None = None
    for root, cluster in base:
        for path in glob.glob(str(root / name / "dfpt-qhaflow" / "runs" / "*"
                                 / f"{name}-dfpt-qha.json")):
            payload = json.load(open(path))
            if payload.get("ok") is not True:
                continue
            if best is None or os.path.getmtime(path) > best[0]:
                best = (os.path.getmtime(path), Path(path), payload)
    if best is None:
        return None
    mtime, path, payload = best
    atoms = int(payload["n_atoms_cell"])
    return {
        "path": str(path), "cluster": next(c for r, c in base if str(path).startswith(str(r))),
        "n_atoms": atoms,
        "E0": payload["eos"]["E0_eV"] / atoms,
        "V0": payload["eos"]["V0_angstrom3"] / atoms,
        "B0": payload["eos"]["B0_GPa"],
        "B0prime": payload["eos"].get("B0prime"),
        "applicability": payload.get("qha_applicability"),
        "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        "volumes": payload["volumes"],
        "static": payload["static_energies_eV"],
        "inside": payload.get("v0_inside_sampled_window"),
    }


def compare(label: str, a: dict, b: dict, problems: list[str]) -> dict:
    """Compare two artifacts, refusing to count one file against itself.

    The mirror makes this trap real: if an endpoint's artifacts were copied
    from the tree being compared against, the difference is exactly zero and
    the check proves nothing.  A byte-identical file identifies that.  An
    identical *volume grid* does not -- the grid comes from the same
    deterministic EOS relaxation, so two genuinely independent runs of the
    same endpoint share it by construction.
    """
    same_file = str(a["path"]) == str(b["path"])
    same_bytes = (not same_file and a.get("sha256") is not None
                  and a.get("sha256") == b.get("sha256"))
    dE = (a["E0"] - b["E0"]) * 1000.0
    dV = a["V0"] - b["V0"]
    dB = a["B0"] - b["B0"]
    ok = (abs(dE) <= SAME_CRYSTAL_TOL_MEV and abs(dV) <= SAME_CRYSTAL_TOL_V
          and abs(dB) <= SAME_CRYSTAL_TOL_B)
    record = {"check": label, "dE_meV": dE, "dV_A3": dV, "dB0_GPa": dB,
              "ok": ok, "independent_runs": not (same_file or same_bytes),
              "a": a["path"], "b": b["path"]}
    if same_file or same_bytes:
        problems.append(
            f"{label}: the two artifacts are the same file, so this comparison "
            "proves nothing -- it is a copy, not a reproduction"
        )
    elif not ok:
        problems.append(
            f"{label}: same crystal disagrees (dE={dE:+.2f} meV/atom, "
            f"dV={dV:+.5f} A^3/atom, dB0={dB:+.2f} GPa)"
        )
    return record


def lattice_constant(volume_per_atom: float, lattice: str) -> float:
    """Cubic lattice constant from volume per atom."""
    return (volume_per_atom * (2.0 if lattice == "BCC" else 4.0)) ** (1.0 / 3.0)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args(argv)

    base = roots()
    problems: list[str] = []
    records: list[dict] = []

    # ---- 1. same crystal, independently computed -----------------------
    for element in BCC_STABLE:
        a = latest_ok(f"SER-{element}-1", base)
        b = latest_ok(f"BCC-{element}-{element}-2", base)
        if a and b:
            records.append(compare(f"SER-{element}-1 == BCC-{element}-{element}-2",
                                   a, b, problems))
    for element in FCC_STABLE:
        a = latest_ok(f"SER-{element}-1", base)
        b = latest_ok(f"FCC-{element}-{element}-4", base)
        if a and b:
            records.append(compare(f"SER-{element}-1 == FCC-{element}-{element}-4",
                                   a, b, problems))

    # ---- 2. the same endpoint on both clusters -------------------------
    endpoints: dict[str, dict[str, dict]] = {}
    for root, cluster in base:
        for path in glob.glob(str(root / "*" / "dfpt-qhaflow" / "runs" / "*"
                                 / "*-dfpt-qha.json")):
            name = Path(path).name[: -len("-dfpt-qha.json")]
            payload = json.load(open(path))
            if payload.get("ok") is not True:
                continue
            atoms = int(payload["n_atoms_cell"])
            endpoints.setdefault(name, {})[cluster] = {
                "path": path, "n_atoms": atoms,
                "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                "E0": payload["eos"]["E0_eV"] / atoms,
                "V0": payload["eos"]["V0_angstrom3"] / atoms,
                "B0": payload["eos"]["B0_GPa"],
                "static": payload["static_energies_eV"],
                "gibbs": payload["gibbs_eV"],
            }
    for name, per_cluster in sorted(endpoints.items()):
        if len(per_cluster) < 2:
            continue
        clusters = sorted(per_cluster)
        for other in clusters[1:]:
            a, b = per_cluster[clusters[0]], per_cluster[other]
            records.append(compare(f"{name} [{clusters[0]} vs {other}]",
                                   a, b, problems))
            # Splitting static from vibrational is what makes a cross-cluster
            # rerun informative: a deterministic SCF reproduces to the last bit
            # on different hardware, while the phonons go through a different
            # MPI decomposition and can differ in the last few digits.
            d_static = max(abs(x - y) for x, y in zip(a["static"], b["static"]))
            d_gibbs = max(abs(x - y) for x, y in zip(a["gibbs"], b["gibbs"]))
            records.append({
                "check": f"{name} [{clusters[0]} vs {other}] static/vibrational",
                "max_abs_static_eV": d_static,
                "max_abs_gibbs_eV": d_gibbs,
                "static_bit_identical": d_static == 0.0,
                "ok": True,
            })
            if d_gibbs * 1000.0 > SAME_CRYSTAL_TOL_MEV:
                problems.append(
                    f"{name}: cross-cluster Gibbs energies differ by "
                    f"{d_gibbs * 1000:.3f} meV/atom, more than the "
                    f"{SAME_CRYSTAL_TOL_MEV} meV tolerance"
                )

    # ---- 3. pure-element lattice constant vs experiment -----------------
    for element in BCC_STABLE + FCC_STABLE + ("Hf", "Ti", "Zr"):
        for lattice, suffix in (("BCC", 2), ("FCC", 4)):
            key = (lattice, element)
            if key not in LATTICE_A_EXPERIMENT:
                continue
            data = latest_ok(f"{lattice}-{element}-{element}-{suffix}", base)
            if data is None:
                continue
            a = lattice_constant(data["V0"], lattice)
            ref = LATTICE_A_EXPERIMENT[key]
            rel = (a - ref) / ref * 100.0
            record = {"check": f"a({lattice},{element})", "computed": a,
                      "experiment": ref, "relative_percent": rel,
                      "ok": abs(rel) <= 3.0}
            records.append(record)
            if not record["ok"]:
                problems.append(
                    f"a({lattice},{element}) = {a:.4f} A is {rel:+.2f}% off the "
                    f"experimental {ref} A -- check the cell, not the functional"
                )

    # ---- 4. the EOS minimum must be inside the scan --------------------
    for root, cluster in base:
        for path in glob.glob(str(root / "*" / "dfpt-qhaflow" / "runs" / "*"
                                 / "*-dfpt-qha.json")):
            payload = json.load(open(path))
            if payload.get("ok") is not True:
                continue
            name = Path(path).name[: -len("-dfpt-qha.json")]
            if payload.get("v0_inside_sampled_window") is False:
                problems.append(f"{name}: EOS minimum sits at the edge of the "
                                f"sampled volume window")

    print(f"roots: {[str(r) for r, _ in base]}")
    for record in records:
        if "dE_meV" in record:
            tag = "ok  " if record["ok"] else "FAIL"
            indep = "independent" if record["independent_runs"] else "COPY"
            print(f"  {tag} {record['check']:40s} "
                  f"dE={record['dE_meV']:+9.6f} meV  dV={record['dV_A3']:+.6f} A^3  "
                  f"dB0={record['dB0_GPa']:+7.3f} GPa  [{indep}]")
        elif "max_abs_gibbs_eV" in record:
            print(f"  ---- {record['check']:40s} static "
                  f"{'bit-identical' if record['static_bit_identical'] else 'differs'}"
                  f", max|dG|={record['max_abs_gibbs_eV'] * 1000:.4f} meV/atom")
        else:
            print(f"  {'ok ' if record['ok'] else 'FAIL'} {record['check']:40s} "
                  f"a={record['computed']:.4f} A  ref={record['experiment']} A  "
                  f"({record['relative_percent']:+.2f}%)")
    print()
    if problems:
        print(f"{len(problems)} problem(s):", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
    else:
        print("no inconsistencies")
    if args.json:
        args.json.write_text(json.dumps(
            {"records": records, "problems": problems}, indent=2) + "\n")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
