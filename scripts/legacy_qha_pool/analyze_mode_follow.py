#!/usr/bin/env python
"""Analyse the minimal mode-following batch (report 4.13.11).

Three things per case, and the ORDER matters:

1. ELECTRONIC + GEOMETRIC GATE, then
2. ENERGY AND MODE-PROJECTED FORCE, then
3. (only for paths that actually descend) the relaxation.

Step 2 deliberately does NOT stop at the total energy.  The total energy of a
single point tells you almost nothing on its own -- a 0.02 A distortion changes
it by ~0.2 meV/atom, the same order as the k-mesh noise.  What identifies a real
soft mode is the sign of the curvature AND the force pushing the structure away
from the parent phase.

Sign convention for the mode-projected force
--------------------------------------------
With the displacement field `u_i` (cartesian, Angstrom) and the forces `F_i`:

    u_hat_i = u_i / |u|          (|u| = the whole supercell's displacement norm)
    F_mode  = sum_i F_i . u_hat_i

For a NEGATIVE-curvature mode (what we are looking for):

    dE(d) < 0        and        F_mode * d > 0

i.e. the force pushes further along the direction we already moved, which is the
same statement as negative curvature written in forces.  Both are reported, and
`consistent_with_negative_curvature` is the conjunction, so a sign slip in either
one shows up as a disagreement rather than as a plausible-looking number.

`F_mode` is reported per unit displacement norm; `F_mode * d` is the quantity
with the sign to read.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

from analyze_align import parse_outcar
from analyze_stage0b import read_forces

HARTREE_eV = 27.211386245988
BOHR_ANG = 0.529177210903


def mode_projection(forces: np.ndarray, u: np.ndarray) -> dict:
    """Project the forces onto the applied displacement field.

    `forces` and `u` are both (n_atoms, 3), forces in eV/Angstrom and u in
    Angstrom.  Returns F_mode (eV/Angstrom per unit displacement norm) plus the
    quantities needed to check the sign by hand.
    """
    un = float(np.linalg.norm(u))
    if un < 1e-12:
        return {"error": "the displacement field is zero"}
    u_hat = u / un
    # F . u_hat summed over atoms: units eV/Angstrom
    f_mode = float((forces * u_hat).sum())
    return {
        "u_norm_angstrom": un,
        "max_single_atom_displacement_angstrom": float(
            np.linalg.norm(u, axis=1).max()),
        "F_mode_eV_per_angstrom": f_mode,
        "F_mode_times_displacement": f_mode * un,
        "F_mode_definition": ("sum_i F_i . u_hat_i with u_hat = u/|u|, u in "
                              "Angstrom, F in eV/Angstrom"),
        "per_atom_force_dot_u_hat": [
            float(np.dot(forces[i], u_hat[i])) for i in range(len(forces))],
    }


def last_oszicar_energy(case: Path) -> float | None:
    """E0 from the LAST converged line of OSZICAR.

    Cross-check for the energy column.  `energy(sigma->0)` appears once per SCF
    iteration AND once in the final block; taking the FIRST occurrence instead
    of the last silently returns iteration 1's energy (a positive number near
    the initial-charge-density guess) rather than the converged one.  That is
    the same failure as report 4.6/4.10 -- wrong energy column -- one level
    down: right column, wrong line.  This function gives an independent value
    to compare against, so the mistake shows up as a mismatch.
    """
    osz = case / "OSZICAR"
    if not osz.is_file():
        return None
    hits = re.findall(r"E0=\s*(-?\d+\.\d+)", osz.read_text(errors="replace"))
    return float(hits[-1]) if hits else None


def analyse_case(case: Path, ref_energy: float | None) -> dict:
    meta = json.loads((case / "case.json").read_text())
    out = parse_outcar(case / "OUTCAR")
    f = read_forces(case / "OUTCAR")
    osz_e = last_oszicar_energy(case)
    u = np.array(meta["displacement_vectors_cartesian"], dtype=float)
    rec = {
        "case": case.name,
        "pattern": meta["pattern"],
        "amp": meta["amplitude_target_angstrom"],
        "dir": str(case),
        "finished": out["finished"],
        "scf_iterations": out["n_scf_iterations"],
        "energy_sigma0_eV": out["energy_sigma0_eV"],
        "energy_from_last_oszicar_e0_eV": osz_e,
        "energy_columns_agree": (
            None if (osz_e is None or out["energy_sigma0_eV"] is None)
            else abs(osz_e - out["energy_sigma0_eV"]) < 1e-5),
        "total_magnetization_muB": out["total_magnetization_muB"],
        "max_force_eV_per_angstrom": (float(np.abs(f).max())
                                      if f is not None else None),
        "n_atoms": int(u.shape[0]),
        "max_single_atom_displacement_angstrom":
            meta["max_single_atom_displacement_angstrom"],
    }
    if f is not None:
        rec.update(mode_projection(f, u))
    if ref_energy is not None and out["energy_sigma0_eV"] is not None:
        d = out["energy_sigma0_eV"] - ref_energy
        rec["delta_E_eV"] = d
        rec["delta_E_meV_per_atom"] = 1e3 * d / rec["n_atoms"]
    return rec


def fit_curvature(rows: list[dict]) -> dict:
    """Quadratic coefficient of dE(d) for one pattern.

    dE = 0.5 * k * d^2 + O(d^3); k < 0 is the soft-mode signature.  Reported
    with the residuals so a fit that is not actually quadratic is visible --
    a straight line through three points is not evidence of a parabola.
    """
    xs = np.array([r["max_single_atom_displacement_angstrom"] for r in rows])
    ys = np.array([r["delta_E_meV_per_atom"] for r in rows])
    if len(xs) < 3:
        return {"error": "need at least 3 amplitudes to fit a curvature"}
    A = np.vstack([0.5 * xs ** 2, xs, np.ones_like(xs)]).T
    coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
    pred = A @ coef
    return {
        # coef[0] is the coefficient of 0.5*d^2, so it IS the curvature k in
        # dE = 0.5*k*d^2.  Do NOT double it.
        "k_meV_per_atom_per_angstrom2": float(coef[0]),
        "linear_term": float(coef[1]),
        "offset": float(coef[2]),
        "residuals_meV_per_atom": [float(y - p) for y, p in zip(ys, pred)],
        "max_abs_residual_meV_per_atom": float(np.abs(ys - pred).max()),
        "is_negative_curvature": bool(coef[0] < 0),
        "amplitudes_angstrom": [float(x) for x in xs],
        "delta_E_meV_per_atom": [float(y) for y in ys],
        "note": ("three amplitudes fit three parameters exactly, so the "
                 "residuals are NOT a goodness-of-fit test -- they only catch "
                 "a grossly non-quadratic shape.  Treat k as an estimate."),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    # Accept either the batch root or the mode_follow directory itself, so an
    # archived copy can be analysed without renaming anything.
    root = args.root / "mode_follow"
    if not root.is_dir():
        root = args.root
    if not (root / "validation.json").is_file():
        raise SystemExit(f"no mode_follow batch at {root}")

    rows: list[dict] = []
    skipped = []
    for case in sorted(root.glob("*/*/")):
        if any(p.startswith("_") for p in case.relative_to(root).parts):
            continue
        if not (case / "OUTCAR").is_file():
            continue
        r = analyse_case(case, None)
        if r["energy_sigma0_eV"] is None:
            # OUTCAR exists but has no final energy block yet -> still running
            skipped.append(str(case.relative_to(root)))
            continue
        rows.append(r)
    if skipped:
        print(f"跳过仍在运行/无最终能量的算例 {len(skipped)} 个: {skipped}")

    # E(0) comes from `reference/amp0`, generated with the SAME ISYM=0 /
    # KPOINTS / smearing / POTCAR.  It must NOT be borrowed from a Stage 0B
    # static: those used ISYM=2 and a different irreducible k-point set (75 vs
    # 260), which would put a k-sampling difference straight into the dE we are
    # trying to resolve to ~0.2 meV/atom.
    ref_row = next((r for r in rows if r["pattern"] == "reference"), None)
    ref = ref_row["energy_sigma0_eV"] if ref_row else None
    if ref is None:
        print("!! 没有 reference/amp0 —— dE 无法计算。")
        print("   E(0) 必须来自同 ISYM=0 / 同 KPOINTS 的零位移算例，")
        print("   不能借用 Stage 0B 的 ISYM=2 静态能量。")
    for r in rows:
        if r["pattern"] == "reference":
            r["delta_E_eV"] = 0.0
            r["delta_E_meV_per_atom"] = 0.0
            continue
        if ref is not None:
            r["delta_E_eV"] = r["energy_sigma0_eV"] - ref
            r["delta_E_meV_per_atom"] = 1e3 * r["delta_E_eV"] / r["n_atoms"]
            if "F_mode_times_displacement" in r:
                # The two independent readings of the same curvature:
                #   from the total energy:  dE
                #   from the forces:        -(1/2) sum_i F_i . u_i
                # For a harmonic mode these must agree.  They use DIFFERENT
                # parts of the OUTCAR (energy block vs force block) and the
                # displacement field, so agreement is a real cross-check that
                # the energy column, the forces and the field are all the ones
                # we think they are.
                half_fu = -0.5 * r["F_mode_times_displacement"]
                r["delta_E_from_forces_eV"] = half_fu
                r["delta_E_energy_vs_forces_ratio"] = (
                    r["delta_E_eV"] / half_fu if abs(half_fu) > 1e-12 else None)
                r["consistent_with_negative_curvature"] = bool(
                    r["delta_E_eV"] < 0
                    and r["F_mode_times_displacement"] > 0)

    print(f"{'pattern':>26} {'amp/A':>6} {'E(s0)/eV':>14} {'dE meV/at':>11} "
          f"{'F_mode':>10} {'F_mode*d':>10} {'dE/(-½F·u)':>11} {'neg curv?':>10}")
    print("-" * 108)
    for r in sorted(rows, key=lambda x: (x["pattern"], x["amp"])):
        if r["pattern"] == "reference":
            print(f"{'reference':>26} {0.0:>6.2f} {r['energy_sigma0_eV']:>14.8f} "
                  f"{0.0:>11.4f} {'-':>10} {'-':>10} {'-':>10}")
            continue
        rr = r.get("delta_E_energy_vs_forces_ratio")
        print(f"{r['pattern']:>26} {r['amp']:>6.2f} "
              f"{r['energy_sigma0_eV']:>14.8f} "
              f"{r.get('delta_E_meV_per_atom', float('nan')):>11.4f} "
              f"{r.get('F_mode_eV_per_angstrom', float('nan')):>10.4f} "
              f"{r.get('F_mode_times_displacement', float('nan')):>10.4f} "
              f"{('%.4f' % rr) if rr is not None else '-':>11} "
              f"{str(r.get('consistent_with_negative_curvature', '-')):>10}")

    bad_e = [r for r in rows if r.get("energy_columns_agree") is False]
    if bad_e:
        print("\n!! 能量列交叉检查失败（OUTCAR 末值与 OSZICAR 末个 E0 不一致）：")
        for r in bad_e:
            print(f"   {r['case']}: outcar={r['energy_sigma0_eV']} "
                  f"oszicar={r['energy_from_last_oszicar_e0_eV']}")
    else:
        print("\n能量列交叉检查: 全部通过（OUTCAR 末值 == OSZICAR 末个 E0）")

    if ref_row is not None:
        print(f"\n参考态残余力: max|F| = "
              f"{ref_row.get('max_force_eV_per_angstrom')} eV/Å "
              f"(应接近 0；明显非零说明未位移的胞本身不平衡)")
        print(f"参考态 SCF 步数 = {ref_row.get('scf_iterations')}, "
              f"finished = {ref_row.get('finished')}")

    print()
    curv = {}
    for pat in sorted({r["pattern"] for r in rows}):
        sub = [r for r in rows
               if r["pattern"] == pat and r["pattern"] != "reference"]
        if not sub:
            continue
        if all("delta_E_meV_per_atom" in r for r in sub):
            curv[pat] = fit_curvature(sub)
            c = curv[pat]
            if "error" in c:
                print(f"{pat}: 曲率无法拟合（{c['error']}）——"
                      f"当前只有 {len(sub)} 个振幅，等 0.10 Å 跑完再判")
                continue
            # A/B: negative curvature AND consistent outward force at the SMALL
            # amplitudes.  The small-amplitude behaviour is the discriminator;
            # 0.10 A can descend for reasons that have nothing to do with the
            # harmonic mode.
            small = [r for r in sub if r["amp"] <= 0.05
                     and "consistent_with_negative_curvature" in r]
            force_ok = bool(small) and all(
                r["consistent_with_negative_curvature"] for r in small)
            if c["is_negative_curvature"] and force_ok:
                verdict = "A: 小振幅负曲率 + 力投影一致 -> 可进 ISIF=2 弛豫"
            elif c["is_negative_curvature"] and not force_ok:
                verdict = ("B: 能量下降但力投影不一致 -> **先停**，"
                           "查模式向量/原子顺序/F_mode 符号/u_hat 归一化/"
                           "E(0) 是否 ISYM=0 匹配/能量列")
            else:
                verdict = "C: 无负曲率 -> 该图案不产生畸变相"
            curv[pat]["small_amplitude_force_consistent"] = force_ok
            curv[pat]["verdict"] = verdict
            print(f"{pat}: k = {c['k_meV_per_atom_per_angstrom2']:+.3f} "
                  f"meV/atom/A^2  负曲率: {c['is_negative_curvature']}")
            print(f"    {verdict}")
        else:
            print(f"{pat}: 缺少参考能 E(0)，无法拟合曲率")
    print("\n注意：三个振幅拟合三个参数是精确拟合，残差不是拟合优度检验。")

    if args.json:
        args.json.write_text(json.dumps(
            {"cases": rows, "curvature": curv}, indent=2) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
