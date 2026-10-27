#!/usr/bin/env python
"""E1 analysis: compare DFPT and finite-displacement phonons for BCC-Al-Ta-2.

Both routes were run on the same 16-atom B2 supercell with the same k-mesh,
ENCUT, smearing and EDIFF. DFPT reduced to 2 perturbations (atoms 1 and 9,
direction x) and phonopy reduced to 2 displacements. So the two datasets carry
the same amount of information and can be fed through the same machinery.

The only thing that differs is how the force response was obtained:
  DFPT  - VASP solved the Sternheimer equations; DYNMAT holds the response to a
          UNIT displacement (so the numbers are force constants already).
  FD    - phonopy displaced by 0.01 A and VASP returned forces; dividing by
          the displacement gives the force constants.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

ROOT = Path("/nfs_ssd/tmp/legacy-qha-pool")
DFPT = ROOT / "e1-dfpt-al-ta"
FD = ROOT / "e1-fdisp-al-ta"
A0 = 3.2246015525754168
SUPER = 2
N = 16


def parse_poscar(p: Path):
    lines = p.read_text().splitlines()
    scale = float(lines[1])
    lat = np.array([[float(x) for x in lines[2 + i].split()[:3]] for i in range(3)]) * scale
    names = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    start = 8 if lines[7].strip().lower().startswith("d") else 8
    frac = np.array([[float(x) for x in lines[start + i].split()[:3]] for i in range(sum(counts))])
    syms = []
    for nm, c in zip(names, counts):
        syms += [nm] * c
    return lat, frac, syms


def parse_dynmat(p: Path):
    """-> list of (atom_index_1based, direction, forces[N,3]) in eV/A^2 units."""
    lines = p.read_text().splitlines()
    n_pert = int(lines[0].split()[-1])
    out = []
    i = 2
    for _ in range(n_pert):
        head = lines[i].split()
        atom, direction = int(head[0]), int(head[1])
        i += 1
        f = np.array([[float(x) for x in lines[i + k].split()[:3]] for k in range(N)])
        i += N
        out.append({"atom": atom, "direction": direction, "forces": f})
    return out


def parse_outcar_forces(p: Path) -> np.ndarray:
    txt = p.read_text(errors="replace")
    # last TOTAL-FORCE block
    blocks = re.findall(
        r"POSITION\s+TOTAL-FORCE \(eV/Angst\)\s*-+\n(.*?)\n\s*-{5,}", txt, re.S)
    if not blocks:
        raise SystemExit(f"no TOTAL-FORCE block in {p}")
    rows = [l.split() for l in blocks[-1].strip().splitlines()]
    return np.array([[float(x) for x in r[3:6]] for r in rows])


def fc_from_forces(cell: PhonopyAtoms, distances: list[float],
                   forces_per_disp: list[np.ndarray], disp_meta: list[dict],
                   supercell_matrix) -> Phonopy:
    """Build a Phonopy object by replaying displacements and their forces.

    Phonopy's own generate_displacements() is called first so that the
    displacement bookkeeping (which atom, which direction, symmetry expansion)
    is exactly the one that produced the POSCARs, then the forces are attached.
    """
    ph = Phonopy(cell, supercell_matrix, primitive_matrix="auto", is_symmetry=True)
    ph.generate_displacements(distance=distances[0])
    ph.forces = forces_per_disp
    ph.produce_force_constants()
    return ph


def main():
    lat_u, frac_u, syms_u = parse_poscar(FD / "POSCAR-unitcell")
    cell = PhonopyAtoms(symbols=syms_u, cell=lat_u, scaled_positions=frac_u)
    scm = np.diag([SUPER, SUPER, SUPER])

    # --- finite displacement -------------------------------------------------
    _, frac_fd, syms_fd = parse_poscar(FD / "disp-001" / "POSCAR")
    fd_forces = [parse_outcar_forces(FD / f"disp-{i:03d}" / "OUTCAR") for i in (1, 2)]
    ph_fd = Phonopy(cell, scm, primitive_matrix="auto", is_symmetry=True)
    ph_fd.generate_displacements(distance=0.01)
    assert len(ph_fd.supercells_with_displacements) == 2
    ph_fd.forces = fd_forces
    ph_fd.produce_force_constants()
    fc_fd = ph_fd.force_constants
    d_fd = ph_fd.dataset

    # --- DFPT ----------------------------------------------------------------
    dyn = parse_dynmat(DFPT / "DYNMAT")
    # VASP's atom order is my POSCAR's order (species-major); phonopy's supercell
    # order is whatever it built. Map by position.
    _, frac_vasp, syms_v = parse_poscar(DFPT / "POSCAR")
    ph_ref = Phonopy(cell, scm, primitive_matrix="auto", is_symmetry=True)
    frac_ph = ph_ref.supercell.scaled_positions
    syms_ph = ph_ref.supercell.symbols
    order = []
    for k in range(N):
        d = np.abs((frac_vasp - frac_ph[k] + 0.5) % 1.0 - 0.5).sum(axis=1)
        j = int(np.argmin(d))
        assert d[j] < 1e-6, (k, d[j])
        assert syms_v[j] == syms_ph[k]
        order.append(j)

    ph_d = Phonopy(cell, scm, primitive_matrix="auto", is_symmetry=True)
    ph_d.generate_displacements(distance=1.0)   # unit displacement: DYNMAT is already Phi
    assert len(ph_d.supercells_with_displacements) == 2
    # reorder DFPT forces into phonopy's supercell order
    dfpt_forces = [d["forces"][order] for d in dyn]
    ph_d.forces = dfpt_forces
    ph_d.produce_force_constants()
    fc_d = ph_d.force_constants

    # --- compare the force constants ----------------------------------------
    # both are (n_super, n_super, 3, 3) with the same atom order now
    diff = fc_d - fc_fd
    scale = np.abs(fc_fd).max()
    print("=== FORCE CONSTANTS (eV/A^2) ===")
    print(f"  max |Phi|        DFPT {np.abs(fc_d).max():10.4f}   FD {np.abs(fc_fd).max():10.4f}")
    print(f"  max |dPhi|       {np.abs(diff).max():10.4f}   "
          f"({100 * np.abs(diff).max() / scale:.2f}% of max |Phi|)")
    print(f"  ASR residual     DFPT {np.abs(fc_d.sum(axis=1)).max():10.4f}   "
          f"FD {np.abs(fc_fd.sum(axis=1)).max():10.4f}")
    print(f"  permutation dev  DFPT {np.abs(fc_d - fc_d.transpose(1,0,3,2)).max():10.4f}   "
          f"FD {np.abs(fc_fd - fc_fd.transpose(1,0,3,2)).max():10.4f}")

    # per-atom-pair comparison for the row of the displaced Al atom
    print("\n  on-site and neighbour Phi_xx of the displaced Al atom (eV/A^2):")
    print(f"    {'pair':>5} {'DFPT':>12} {'FD':>12} {'ratio':>8}")
    for j in range(N):
        a, b = fc_d[0, j, 0, 0], fc_fd[0, j, 0, 0]
        r = a / b if abs(b) > 1e-12 else float("nan")
        print(f"    {j:>5} {a:>12.5f} {b:>12.5f} {r:>8.3f}")

    # --- compare the phonons -------------------------------------------------
    print("\n=== PHONONS ===")
    out = {"schema": "legacy-qha-pool/e1-comparison/v1",
           "endpoint": "BCC-Al-Ta-2", "cell": "B2 AlTa 2x2x2, a0=3.2246015525754168",
           "force_constants": {
               "dfpt_max_abs": float(np.abs(fc_d).max()),
               "fd_max_abs": float(np.abs(fc_fd).max()),
               "max_abs_diff": float(np.abs(diff).max()),
               "max_abs_diff_percent_of_scale": float(100 * np.abs(diff).max() / scale),
               "dfpt_asr_residual": float(np.abs(fc_d.sum(axis=1)).max()),
               "fd_asr_residual": float(np.abs(fc_fd.sum(axis=1)).max()),
               "dfpt_permutation_dev": float(np.abs(fc_d - fc_d.transpose(1,0,3,2)).max()),
               "fd_permutation_dev": float(np.abs(fc_fd - fc_fd.transpose(1,0,3,2)).max()),
           }, "frequencies": {}}

    for name, ph in (("DFPT", ph_d), ("FD", ph_fd)):
        ph.run_qpoints([[0.0, 0.0, 0.0],
                        [0.5, 0.5, 0.0],
                        [0.5, 0.5, 0.5]])
        fr = ph.get_qpoints_dict()["frequencies"]
        g, Npt, H = (list(map(float, f)) for f in fr)
        print(f"  {name:5s} Gamma {[round(x,4) for x in g]}")
        print(f"        N(1/2,1/2,0) {[round(x,4) for x in Npt]}")
        print(f"        H(1/2,1/2,1/2) {[round(x,4) for x in H]}")
        out["frequencies"][name] = {"Gamma": g, "N": Npt, "H": H}
        ph.run_mesh([12, 12, 12])
        ph.run_thermal_properties(t_min=0, t_max=1000, t_step=100)
        tp = ph.thermal_properties
        print(f"        ZPE {tp.zero_point_energy:.6f} kJ/mol   "
              f"modes {tp.number_of_integrated_modes}/{tp.number_of_modes}")
        out["frequencies"][name]["zpe_kJ_per_mol"] = float(tp.zero_point_energy)
        out["frequencies"][name]["modes_total"] = int(tp.number_of_modes)
        out["frequencies"][name]["modes_integrated"] = int(tp.number_of_integrated_modes)

    dg = np.array(out["frequencies"]["DFPT"]["Gamma"]) - np.array(out["frequencies"]["FD"]["Gamma"])
    out["frequency_delta_THz"] = {"Gamma": [float(x) for x in dg],
                                  "max_abs_GHz": float(np.abs(dg).max() * 1000)}
    print(f"\n  Delta Gamma (DFPT - FD): {[round(x,4) for x in dg]} THz  "
          f"(max {np.abs(dg).max()*1000:.1f} GHz)")
    (ROOT / "e1-comparison.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"\n-> {ROOT/'e1-comparison.json'}")


if __name__ == "__main__":
    main()
