#!/usr/bin/env python
"""Stage 0B analysis: forces -> force constants -> phonons, across smearings.

Three levels, in the order the review asked for:

  1. FORCES        per displacement: max|F|, RMS, per-atom force differences
                   between smearings, SCF iterations, magnetic state.
  2. FORCE         build FORCE_CONSTANTS with phonopy from each smearing's
     CONSTANTS     displacements, then compare max|dPhi|, relative difference,
                   ASR residual (absolute and relative) and the permutation
                   residual.  A single matrix element is not enough -- it can
                   agree by accident while the dispersion does not.
  3. PHONONS       Gamma acoustic/optical, N = (1/2,1/2,0), H = (1/2,1/2,1/2),
                   X = (0,1/2,0), a [110] point, the minimum frequency on a
                   mesh, the count of q-points with imaginary frequencies, ZPE.

Ordering safety
---------------
The force array handed to phonopy must be in phonopy's own supercell order.  The
POSCARs were written from phonopy's supercell (build_stage0b.phonopy_supercell),
and every case records a ``cell_fingerprint``.  This script rebuilds the same
phonopy object, asserts the fingerprint matches, and asserts that each POSCAR
differs from the undisplaced supercell by exactly the recorded displacement on
exactly the recorded atom.  Without those two asserts a permutation would go
unnoticed: the force constants would still be nonzero, roughly symmetric and
give plausible frequencies.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.core import Structure

import build_stage0 as S0
import build_stage0b as B0
from analyze_align import parse_outcar
from analyze_stage0 import parse_site_moments

SMEARINGS = ("s1_010", "s1_020", "s2_010", "s2_020")
SMEARINGS_MAP = {"s1_010": (1, 0.10), "s1_020": (1, 0.20),
                 "s2_010": (2, 0.10), "s2_020": (2, 0.20)}
REF = "s1_020"
LABEL = {"s1_010": "1/0.10", "s1_020": "1/0.20", "s2_010": "2/0.10", "s2_020": "2/0.20"}
AMPS = ("0.01", "0.03")
QPOINTS = {"Gamma": [0.0, 0.0, 0.0], "X": [0.0, 0.5, 0.0],
           "N": [0.5, 0.5, 0.0], "N110_0.4": [0.4, 0.4, 0.0],
           "H": [0.5, 0.5, 0.5]}


def smearing_tag(meta: dict) -> str:
    """Recover the smearing tag from a case's recorded ISMEAR/SIGMA.

    The directory name is ``a{amp}_d{disp:02d}_{tag}`` and ``tag`` itself
    contains an underscore (``s1_010``), so ``name.split("_", 1)[1]`` yields
    ``d01_s1_010`` -- which matches nothing and silently produced empty tables.
    Deriving the tag from the recorded ISMEAR/SIGMA instead is unambiguous.
    """
    ism = int(meta["ismear"])
    sig = meta.get("sigma")
    sig = None if sig in (None, "None") else float(sig)
    for tag, (i, s) in SMEARINGS_MAP.items():
        if i == ism and ((s is None and sig is None) or
                         (s is not None and sig is not None and abs(s - sig) < 1e-9)):
            return tag
    raise KeyError(f"no smearing tag for ISMEAR={ism} SIGMA={sig}")


def read_forces(outcar: Path) -> np.ndarray | None:
    """Per-atom forces (eV/Angstrom) from the LAST TOTAL-FORCE block."""
    txt = outcar.read_text(errors="replace")
    blocks = txt.split("TOTAL-FORCE (eV/Angst)")
    if len(blocks) < 2:
        return None
    rows = []
    for ln in blocks[-1].splitlines():
        parts = ln.split()
        if len(parts) == 6:
            try:
                rows.append([float(x) for x in parts[3:6]])
            except ValueError:
                continue
        elif rows:
            break
    return np.array(rows) if rows else None


def make_phonopy(endpoint: str) -> tuple[Phonopy, Structure]:
    """Same construction build_stage0b used, so the ordering is identical."""
    base = S0.build_structures()[endpoint]["struct"]
    return B0.phonopy_supercell(base, 2)


def build_fc(endpoint: str, dirs: list[Path], out_dir: Path) -> dict:
    """FORCE_CONSTANTS from one smearing's displaced directories.

    phonopy's force setter needs a dataset, so the displacements are regenerated
    with the same amplitude and settings the builder used, and then checked
    against what each case recorded.  That check is the point: if the builder's
    displacement list and this one ever disagree (different phonopy version,
    different symmetry tolerance), the forces would be attached to the wrong
    displacements and the force constants would be silently wrong.
    """
    ph, sup = make_phonopy(endpoint)
    metas = [json.loads((d / "case.json").read_text()) for d in dirs]
    amp = float(metas[0]["amplitude_angstrom"])
    fp = B0.cell_fingerprint(sup)
    for d, meta in zip(dirs, metas):
        if meta.get("cell_fingerprint") != fp:
            return {"error": f"{d}: cell fingerprint mismatch"}
        if abs(float(meta["amplitude_angstrom"]) - amp) > 1e-12:
            return {"error": f"{d}: amplitude differs within one set"}

    # regenerate the displacement set exactly as the builder did
    ph.generate_displacements(distance=amp, is_plusminus="auto", is_diagonal=True)
    if len(ph.displacements) != len(dirs):
        return {"error": f"{endpoint}: regenerated {len(ph.displacements)} "
                         f"displacements, but {len(dirs)} directories were given"}
    for i, (d, meta) in enumerate(zip(dirs, metas)):
        dd = ph.displacements[i]
        vec = [float(dd[1]), float(dd[2]), float(dd[3])]
        if int(dd[0]) != int(meta["disp_atom"]) or \
                not np.allclose(vec, meta["disp_vector"], atol=1e-9):
            return {"error": f"{d}: regenerated displacement {i} is atom {int(dd[0])} "
                             f"{vec}, but the case recorded atom {meta['disp_atom']} "
                             f"{meta['disp_vector']}"}

    # each POSCAR must be `sup` displaced on exactly the recorded atom
    forces = []
    for d, meta in zip(dirs, metas):
        st = Structure.from_file(d / "POSCAR")
        if len(st) != len(sup):
            return {"error": f"{d}: {len(st)} atoms, expected {len(sup)}"}
        delta = st.cart_coords - sup.cart_coords
        moved = np.where(np.linalg.norm(delta, axis=1) > 1e-8)[0]
        if len(moved) != 1 or int(moved[0]) != int(meta["disp_atom"]):
            return {"error": f"{d}: moved atoms {moved.tolist()}, recorded "
                             f"{meta['disp_atom']}"}
        if not np.allclose(delta[moved[0]], meta["disp_vector"], atol=1e-6):
            return {"error": f"{d}: displacement vector does not match the record"}
        f = read_forces(d / "OUTCAR")
        if f is None:
            return {"error": f"no forces in {d}"}
        if f.shape != (len(sup), 3):
            return {"error": f"{d}: forces shape {f.shape}, expected "
                             f"({len(sup)}, 3)"}
        forces.append(f)

    ph.forces = np.array(forces)
    ph.produce_force_constants(show_drift=False)
    fc = np.array(ph.force_constants, dtype=float)     # (n_satom, n_satom, 3, 3)

    asr = fc.sum(axis=1)
    asr_abs = float(np.abs(asr).max())
    scale = float(np.abs(fc).max())
    # permutation: Phi[i,j,a,b] == Phi[j,i,b,a]; compare on the full square block
    perm = float(np.abs(fc - np.transpose(fc, (1, 0, 3, 2))).max())

    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / "fc.npy", fc)
    return {"max_abs_Phi": scale, "asr_abs": asr_abs,
            "asr_rel": asr_abs / scale if scale else float("nan"),
            "permutation_residual": perm, "shape": list(fc.shape),
            "n_displacements": len(dirs), "fc_file": str(out_dir / "fc.npy")}


def phonons(endpoint: str, fc: np.ndarray) -> dict:
    """Frequencies (THz) at explicit q-points plus a mesh scan."""
    ph, _ = make_phonopy(endpoint)
    ph.force_constants = fc
    res: dict = {}
    ph.run_qpoints([QPOINTS[k] for k in QPOINTS])
    freqs = np.array(ph.get_qpoints_dict()["frequencies"])
    for i, name in enumerate(QPOINTS):
        res[name] = [float(x) for x in freqs[i]]
    ph.run_mesh([12, 12, 12], with_eigenvectors=False)
    mf = np.array(ph.get_mesh_dict()["frequencies"])
    res["mesh_min_THz"] = float(mf.min())
    res["mesh_n_negative"] = int((mf < -1e-4).sum())
    res["mesh_n_modes"] = int(mf.size)
    return res


def zpe(fc: np.ndarray, endpoint: str) -> float | None:
    """Zero-point energy in kJ/mol, from a mesh integration."""
    ph, _ = make_phonopy(endpoint)
    ph.force_constants = fc
    ph.run_mesh([12, 12, 12])
    ph.run_thermal_properties(t_min=0, t_max=100, t_step=100)
    return float(ph.thermal_properties.zero_point_energy)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--skip-fc", action="store_true")
    args = ap.parse_args()
    root = args.root / "stage0b"
    if not root.is_dir():
        raise SystemExit(f"no {root}")

    # ---------------- 1. forces ------------------------------------------
    rows: list[dict] = []
    for case in sorted(root.glob("*/*/")):
        if any(p.startswith(".") for p in case.relative_to(root).parts):
            continue
        if not (case / "OUTCAR").is_file():
            continue
        meta = json.loads((case / "case.json").read_text())
        o = parse_outcar(case / "OUTCAR")
        f = read_forces(case / "OUTCAR")
        rows.append({
            "case": case.name, "endpoint": meta["endpoint"],
            "amp": f"{meta['amplitude_angstrom']:g}",
            "disp": meta["disp_index"], "disp_atom": meta["disp_atom"],
            "disp_species": meta.get("disp_atom_species"),
            "tag": smearing_tag(meta),
            "natoms": meta["n_supercell_atoms"],
            "finished": o["finished"], "scf": o["n_scf_iterations"],
            "energy_sigma0_eV": o["energy_sigma0_eV"],
            "mag": o["total_magnetization_muB"],
            "site_moments": parse_site_moments(case / "OUTCAR"),
            "max_force": (float(np.abs(f).max()) if f is not None else None),
            "rms_force": (float(np.sqrt((f ** 2).sum(axis=1).mean()))
                          if f is not None else None),
            "forces": (f.tolist() if f is not None else None),
            "dir": str(case)})

    print(f"{len(rows)} cases with OUTCAR")
    bad = [r for r in rows if not r["finished"] or r["forces"] is None]
    for r in bad:
        print(f"   !! {r['endpoint']}/{r['case']}")
    summary: dict = {"n_cases": len(rows), "n_bad": len(bad), "rows": rows}

    endpoints = sorted({r["endpoint"] for r in rows})

    # ---------------- 1. forces table ------------------------------------
    print("\n=== 1. forces: max|F| / RMS per case (eV/A) ===")
    hdr = (f"{'endpoint':<15} {'amp':>5} {'d':>2} {'atom':>5} "
           + " ".join(f"{LABEL[t]:>17}" for t in SMEARINGS))
    print(hdr)
    print("-" * len(hdr))
    for ep in endpoints:
        for amp in AMPS:
            disps = sorted({r["disp"] for r in rows
                            if r["endpoint"] == ep and r["amp"] == amp})
            for disp in disps:
                cells, atom, sp = [], "?", "?"
                for t in SMEARINGS:
                    r = next((x for x in rows if x["endpoint"] == ep and x["amp"] == amp
                              and x["disp"] == disp and x["tag"] == t), None)
                    if r is None or r["max_force"] is None:
                        cells.append(f"{'-':>17}")
                        continue
                    atom, sp = r["disp_atom"], r["disp_species"]
                    cells.append(f"{r['max_force']:>8.5f}/{r['rms_force']:<8.5f}")
                print(f"{ep:<15} {amp:>5} {disp:>2} {atom}({sp}): "
                      + " ".join(cells))

    print("\n=== 1b. per-atom force difference vs ISMEAR=1/SIGMA=0.20 (eV/A) ===")
    print(f"{'endpoint':<15} {'amp':>5} {'d':>2} "
          + " ".join(f"{LABEL[t]:>19}" for t in SMEARINGS if t != REF))
    for ep in endpoints:
        for amp in AMPS:
            for disp in sorted({r["disp"] for r in rows
                                if r["endpoint"] == ep and r["amp"] == amp}):
                ref = next((x for x in rows if x["endpoint"] == ep and x["amp"] == amp
                            and x["disp"] == disp and x["tag"] == REF), None)
                if ref is None or ref["forces"] is None:
                    continue
                fr = np.array(ref["forces"])
                cells = []
                for t in SMEARINGS:
                    if t == REF:
                        continue
                    r = next((x for x in rows if x["endpoint"] == ep and x["amp"] == amp
                              and x["disp"] == disp and x["tag"] == t), None)
                    if r is None or r["forces"] is None:
                        cells.append(f"{'-':>19}")
                        continue
                    df = np.abs(np.array(r["forces"]) - fr)
                    cells.append(f"{df.max():>9.5f}/{np.sqrt((df**2).mean()):<9.5f}")
                print(f"{ep:<15} {amp:>5} {disp:>2} " + " ".join(cells))

    # ---------------- 2. force constants ---------------------------------
    fc_summary: dict = {}
    cache: dict = {}
    if not args.skip_fc:
        print("\n=== 2. force constants per smearing ===")
        for ep in endpoints:
            for amp in AMPS:
                print(f"  {ep}  amp={amp}")
                print(f"    {'smearing':>9} {'max|Phi|':>10} {'ASR abs':>11} "
                      f"{'ASR rel':>11} {'perm':>11} {'nDisp':>6}")
                for t in SMEARINGS:
                    rs = sorted((x for x in rows if x["endpoint"] == ep
                                 and x["amp"] == amp and x["tag"] == t
                                 and x["forces"] is not None),
                                key=lambda x: x["disp"])
                    if not rs:
                        continue
                    dirs = [Path(x["dir"]) for x in rs]
                    out = Path(rs[0]["dir"]).parent / f".fc_{amp}_{t}"
                    info = build_fc(ep, dirs, out)
                    if "error" in info:
                        print(f"    {LABEL[t]:>9} ERROR: {info['error']}")
                        continue
                    cache[(ep, amp, t)] = info["fc_file"]
                    fc_summary.setdefault(ep, {}).setdefault(amp, {})[t] = info
                    print(f"    {LABEL[t]:>9} {info['max_abs_Phi']:>10.6f} "
                          f"{info['asr_abs']:>11.3e} {info['asr_rel']:>11.3e} "
                          f"{info['permutation_residual']:>11.3e} "
                          f"{info['n_displacements']:>6}")
        summary["force_constants"] = fc_summary

        print("\n=== 2b. force-constant difference vs ISMEAR=1/SIGMA=0.20 ===")
        print(f"{'endpoint':<15} {'amp':>5} "
              + " ".join(f"{LABEL[t]:>24}" for t in SMEARINGS if t != REF))
        print("   (max|dPhi| eV/A^2 / as a fraction of the reference max|Phi|)")
        for ep in endpoints:
            for amp in AMPS:
                kref = (ep, amp, REF)
                if kref not in cache:
                    continue
                fc_ref = np.load(cache[kref])
                denom = np.abs(fc_ref).max()
                cells = []
                for t in SMEARINGS:
                    if t == REF:
                        continue
                    k = (ep, amp, t)
                    if k not in cache:
                        cells.append(f"{'-':>24}")
                        continue
                    d = np.abs(np.load(cache[k]) - fc_ref).max()
                    cells.append(f"{d:>11.5f} / {d / denom:<11.4%}")
                print(f"{ep:<15} {amp:>5} " + " ".join(cells))

    # ---------------- 3. phonons ------------------------------------------
    if cache:
        print("\n=== 3. phonons (THz) ===")
        ph_summary: dict = {}
        for ep in endpoints:
            for amp in AMPS:
                for t in SMEARINGS:
                    k = (ep, amp, t)
                    if k not in cache:
                        continue
                    fc = np.load(cache[k])
                    r = phonons(ep, fc)
                    z = zpe(fc, ep)
                    r["zpe_kJ_per_mol"] = z
                    ph_summary.setdefault(ep, {}).setdefault(amp, {})[t] = r
        summary["phonons"] = ph_summary

        for ep in endpoints:
            for amp in AMPS:
                if amp not in ph_summary.get(ep, {}):
                    continue
                print(f"\n  {ep}  amp={amp}")
                print(f"    {'quantity':>14} " + " ".join(f"{LABEL[t]:>22}"
                                                          for t in SMEARINGS))
                for qn in list(QPOINTS) + ["mesh_min_THz", "mesh_n_negative",
                                           "zpe_kJ_per_mol"]:
                    cells = []
                    for t in SMEARINGS:
                        r = ph_summary[ep][amp].get(t)
                        if not r or r.get(qn) is None:
                            cells.append(f"{'-':>22}")
                            continue
                        v = r[qn]
                        if isinstance(v, list):
                            cells.append(f"{min(v):>8.4f}..{max(v):<10.4f}")
                        elif isinstance(v, float):
                            cells.append(f"{v:>22.5f}")
                        else:
                            cells.append(f"{v:>22}")
                    print(f"    {qn:>14} " + " ".join(cells))

    # ---------------- 4. magnetic stability -------------------------------
    print("\n=== 4. magnetic state across displacements and smearings ===")
    mag_summary: dict = {}
    for ep in endpoints:
        print(f"  {ep}")
        for t in SMEARINGS:
            rs = [r for r in rows if r["endpoint"] == ep and r["tag"] == t
                  and r["finished"]]
            if not rs:
                continue
            mags = [r["mag"] for r in rs if r["mag"] is not None]
            sites = [r["site_moments"] for r in rs if r["site_moments"]]
            sr = "n/a"
            if sites:
                flat = list(zip(*sites))
                sr = ", ".join(f"{min(x):+.3f}..{max(x):+.3f}" for x in flat)
            mag_summary.setdefault(ep, {})[t] = {
                "mag_range": [min(mags), max(mags)] if mags else None,
                "site_ranges": sr}
            mr = f"{min(mags):+.4f}..{max(mags):+.4f}" if mags else "n/a"
            print(f"    {LABEL[t]:>9}: tot {mr:>20}   sites {sr}")

    # ---------------- 5. verdict ------------------------------------------
    print("\n=== 5. verdict ===")
    if cache and summary.get("phonons"):
        for ep in endpoints:
            for amp in AMPS:
                cells = []
                for t in SMEARINGS:
                    k = (ep, amp, t)
                    if k not in cache:
                        continue
                    fc_ref = np.load(cache[(ep, amp, REF)]) if (ep, amp, REF) in cache else None
                    d = (np.abs(np.load(cache[k]) - fc_ref).max() / np.abs(fc_ref).max()
                         if fc_ref is not None else float("nan"))
                    r = summary["phonons"].get(ep, {}).get(amp, {}).get(t, {})
                    n_neg = r.get("mesh_n_negative")
                    mn = r.get("mesh_min_THz")
                    cells.append(f"{LABEL[t]}: dPhi {d:.2%}  min {mn:.3f}  "
                                 f"nNeg {n_neg}")
                print(f"  {ep:<15} amp={amp}  " + " | ".join(cells))
        print("\n   reading: if dPhi stays small, the minimum frequency and the")
        print("   imaginary-mode count agree across smearings, the endpoint is")
        print("   phonon_recipe_robust; if the soft-mode classification changes,")
        print("   it is phonon_recipe_sensitive.")

    if args.json:
        args.json.write_text(json.dumps(summary, indent=2, default=str) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
