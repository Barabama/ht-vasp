#!/usr/bin/env python3
"""
Generate CoNi(OH)2 SQS structure from Ni(OH)2 3x2x1 orthorhombic supercell.
Replace half of Ni with Co, use mcsqs to find optimal SQS configuration.

Usage:
  .conda/bin/python gen_conioh_sqs.py
Output:
  CoNiOH2-321.vasp
"""

import os
import sys

# Ensure AT-AT binaries are in PATH
os.environ["PATH"] = "/workspace/softwares/atat/bin:" + os.environ.get("PATH", "")

from pymatgen.core import Structure, Lattice
from pymatgen.command_line.mcsqs_caller import run_mcsqs
from pymatgen.transformations.standard_transformations import SupercellTransformation
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def build_ni_oh2_ortho_supercell(nx=3, ny=2, nz=1):
    """Build Ni(OH)2 hexagonal -> orthorhombic -> nx×ny×nz supercell."""
    a_hex, c_hex = 3.12512265, 4.47344830
    hex_lat = Lattice.hexagonal(a_hex, c_hex)
    species = ["Ni", "H", "H", "O", "O"]
    coords = [
        [0.0, 0.0, 0.0],
        [1 / 3, 2 / 3, 0.55907427],
        [2 / 3, 1 / 3, 0.44092573],
        [1 / 3, 2 / 3, 0.77512013],
        [2 / 3, 1 / 3, 0.22487987],
    ]
    hex_struct = Structure(hex_lat, species, coords)

    # hex -> ortho
    T_ortho = np.array([[1, 0, 0], [1, 2, 0], [0, 0, 1]])
    ortho = SupercellTransformation(T_ortho).apply_transformation(hex_struct)

    # supercell
    T_sc = np.diag([nx, ny, nz])
    return SupercellTransformation(T_sc).apply_transformation(ortho)


def main():
    print("=" * 60)
    print("Generating SQS for CoNi(OH)2 (3x2x1, half Co, half Ni)")
    print("=" * 60)

    # Step 1: Build base 3x2x1 supercell
    print("\n[1] Building Ni(OH)2 3x2x1 orthorhombic supercell...")
    base = build_ni_oh2_ortho_supercell(3, 2, 1)

    # Count Ni sites
    ni_indices = [i for i, site in enumerate(base) if site.specie.symbol == "Ni"]
    n_ni = len(ni_indices)
    n_co = n_ni // 2
    print(f"  Formula: {base.formula}")
    print(f"  Total atoms: {len(base)}")
    print(f"  Ni sites: {n_ni}")
    print(f"  Target: {n_co} Co + {n_ni - n_co} Ni")

    # Step 2: Create disordered structure (partial occupancy on Ni sites)
    print("\n[2] Creating disordered structure with Co0.5Ni0.5 on metal sites...")
    disordered = base.copy()
    for idx in ni_indices:
        disordered[idx] = {"Co": 0.5, "Ni": 0.5}
    print(f"  Is ordered: {disordered.is_ordered}")

    # Save original lattice for post-mcsqs reorientation
    original_lattice = base.lattice

    # Step 3: Run mcsqs
    print("\n[3] Running mcsqs...")
    clusters = {
        2: 6.0,   # pairs within 6 Å
        3: 4.0,   # triplets within 4 Å
    }

    result = run_mcsqs(
        structure=disordered,
        clusters=clusters,
        scaling=1,          # already a supercell
        search_time=1.0,    # 1 minute
        temperature=1.0,
        wr=1.0,
        wn=1.0,
        wd=0.5,
        tol=1e-3,
    )

    sqs_structure = result.bestsqs
    print(f"  Objective function: {result.objective_function}")
    print(f"  Composition: {sqs_structure.composition}")
    print(f"  mcsqs lattice: a={sqs_structure.lattice.a:.3f} b={sqs_structure.lattice.b:.3f} c={sqs_structure.lattice.c:.3f}")

    # Step 4: Reorient to original lattice convention
    # mcsqs permutes axes (str2cif reorientation). Restore original orientation.
    print("\n[4] Restoring original lattice orientation...")
    restored = Structure(
        original_lattice,
        sqs_structure.species,
        sqs_structure.cart_coords,
        coords_are_cartesian=True,
    )
    print(f"  Restored: a={restored.lattice.a:.4f} b={restored.lattice.b:.4f} c={restored.lattice.c:.4f}")

    # Step 5: Sort and write POSCAR
    print("\n[5] Writing CoNiOH2-321.vasp...")
    element_order = ["Co", "Ni", "H", "O"]

    sorted_species = []
    sorted_coords = []
    for elem in element_order:
        for site in restored:
            if site.specie.symbol == elem:
                sorted_species.append(elem)
                sorted_coords.append(site.frac_coords)

    sorted_struct = Structure(restored.lattice, sorted_species, sorted_coords)

    from collections import Counter
    counts = Counter(sorted_species)
    print(f"  Element order: {' '.join(element_order)}")
    print(f"  Counts: {dict(counts)}")

    outpath = os.path.join(SCRIPT_DIR, "CoNiOH2-321.vasp")
    sorted_struct.to(fmt="poscar", filename=outpath)
    print(f"\n  Written to {outpath}")
    print(f"  Lattice: a={sorted_struct.lattice.a:.4f}, b={sorted_struct.lattice.b:.4f}, c={sorted_struct.lattice.c:.4f}")
    print(f"  Angles: {sorted_struct.lattice.angles}")

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
