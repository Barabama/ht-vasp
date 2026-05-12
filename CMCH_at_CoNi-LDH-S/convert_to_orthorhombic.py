#!/usr/bin/env python3
"""
Convert crystal structures to orthorhombic cells for heterojunction modeling.

Supports:
  - Ni(OH)2 (mp-27912): hexagonal P-3m1 → orthorhombic
  - Co2(OH)2CO3 (mp-1202876): monoclinic P2_1/c → orthorhombic → axis reoriented

Usage:
  .conda/bin/python convert_to_orthorhombic.py ni_oh2
  .conda/bin/python convert_to_orthorhombic.py ni_oh2 --nx 3 --ny 2 --nz 2
  .conda/bin/python convert_to_orthorhombic.py co2_oh2_co3
  .conda/bin/python convert_to_orthorhombic.py co2_oh2_co3 --nx 1 --ny 3 --nz 1
"""

import argparse
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import SupercellTransformation


# ---------------------------------------------------------------------------
# Structure builders
# ---------------------------------------------------------------------------

def build_ni_oh2():
    """Ni(OH)2 hexagonal CdI2-type, P-3m1 (#164). mp-27912 params."""
    a, c = 3.12512265, 4.47344830
    lattice = Lattice.hexagonal(a, c)
    species = ["Ni", "H", "H", "O", "O"]
    coords = [
        [0.0, 0.0, 0.0],               # Ni 1a
        [1 / 3, 2 / 3, 0.55907427],    # H 2d
        [2 / 3, 1 / 3, 0.44092573],    # H 2d
        [1 / 3, 2 / 3, 0.77512013],    # O 2d
        [2 / 3, 1 / 3, 0.22487987],    # O 2d
    ]
    return Structure(lattice, species, coords)


def build_co2_oh2_co3():
    """Co2(OH)2CO3 rosasite-type, P2_1/c (#14). mp-1202876 params."""
    a, b, c = 3.19535090, 12.30609000, 9.67280607
    beta = 93.05194829
    lattice = Lattice.from_parameters(a, b, c, 90, beta, 90)

    # Asymmetric unit (Wyckoff 4e, each site ×4 by symmetry)
    asym_species = ["Co", "Co", "H", "H", "C", "O", "O", "O", "O", "O"]
    asym_coords = [
        [0.182247, 0.213507, 0.982709],   # Co0
        [0.490579, 0.112834, 0.714675],   # Co1
        [0.031658, 0.018583, 0.891847],   # H2
        [0.065704, 0.600745, 0.000545],   # H3
        [0.496399, 0.131918, 0.258836],   # C4
        [0.014821, 0.092619, 0.856155],   # O5
        [0.016931, 0.653769, 0.925039],   # O6
        [0.472715, 0.230779, 0.310217],   # O7
        [0.478428, 0.551169, 0.158296],   # O8
        [0.499465, 0.120194, 0.125164],   # O9
    ]

    # P2_1/c symmetry operations (b-unique)
    ops = [
        lambda x, y, z: (x, y, z),
        lambda x, y, z: (-x, -y, -z),
        lambda x, y, z: (-x, y + 0.5, -z + 0.5),
        lambda x, y, z: (x, -y + 0.5, z + 0.5),
    ]

    all_species, all_coords = [], []
    for sp, pos in zip(asym_species, asym_coords):
        x, y, z = pos
        for op in ops:
            fx, fy, fz = op(x, y, z)
            all_species.append(sp)
            all_coords.append([fx % 1.0, fy % 1.0, fz % 1.0])

    return Structure(lattice, all_species, all_coords)


# ---------------------------------------------------------------------------
# Transformation functions
# ---------------------------------------------------------------------------

def hex_to_orthorhombic(structure):
    """Convert hexagonal to orthorhombic via [[1,0,0],[1,2,0],[0,0,1]]."""
    T = np.array([[1, 0, 0], [1, 2, 0], [0, 0, 1]])
    return SupercellTransformation(T).apply_transformation(structure)


def monoclinic_to_orthorhombic(structure):
    """Convert monoclinic to orthorhombic by setting beta=90 (small approx)."""
    a, b, c = structure.lattice.abc
    new_lattice = Lattice.orthorhombic(a, b, c)
    return Structure(
        new_lattice,
        structure.species,
        structure.frac_coords,
        site_properties=structure.site_properties,
    )


def reorient_axes(structure):
    """Reorient: new_a=old_c, new_b=old_a, new_c=old_b.

    Matrix [010, 001, 100] as column vectors, i.e. new = old @ P.
    """
    P = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    return SupercellTransformation(P).apply_transformation(structure)


def make_supercell(structure, nx, ny, nz):
    """Diagonal supercell expansion."""
    T = np.diag([nx, ny, nz])
    return SupercellTransformation(T).apply_transformation(structure)


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def print_info(label, structure):
    """Print structure summary."""
    abc = structure.lattice.abc
    ang = structure.lattice.angles
    sga = SpacegroupAnalyzer(structure)
    print(f"  [{label}]")
    print(f"    Formula:  {structure.formula}")
    print(f"    Atoms:    {len(structure)}")
    print(f"    a,b,c:    {abc[0]:.4f}  {abc[1]:.4f}  {abc[2]:.4f} Å")
    print(f"    α,β,γ:    {ang[0]:.4f}  {ang[1]:.4f}  {ang[2]:.4f}°")
    print(f"    Volume:   {structure.lattice.volume:.3f} Å³")
    print(f"    Space grp: {sga.get_space_group_symbol()} (#{sga.get_space_group_number()})")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert crystal structures to orthorhombic cells"
    )
    parser.add_argument(
        "structure",
        choices=["ni_oh2", "co2_oh2_co3"],
        help="Structure to convert",
    )
    parser.add_argument("--nx", type=int, default=1, help="Supercell multiplier along a")
    parser.add_argument("--ny", type=int, default=1, help="Supercell multiplier along b")
    parser.add_argument("--nz", type=int, default=1, help="Supercell multiplier along c")
    args = parser.parse_args()

    # --- Build and convert ---
    if args.structure == "ni_oh2":
        default_nx, default_ny, default_nz = 3, 2, 2
        # Use defaults if user didn't specify any
        if not any(f"--{ax}" in " ".join(__import__("sys").argv) for ax in ["nx", "ny", "nz"]):
            args.nx, args.ny, args.nz = default_nx, default_ny, default_nz

        print("=" * 60)
        print("Ni(OH)2: hexagonal P-3m1 → orthorhombic")
        print("=" * 60)

        s = build_ni_oh2()
        print_info("1. Hexagonal (original)", s)

        s = hex_to_orthorhombic(s)
        print_info("2. Orthorhombic", s)

        s = make_supercell(s, args.nx, args.ny, args.nz)
        print_info(f"3. Supercell {args.nx}x{args.ny}x{args.nz}", s)

        fname = f"POSCAR_ni_oh2_ortho_{args.nx}x{args.ny}x{args.nz}.vasp"

    else:  # co2_oh2_co3
        default_nx, default_ny, default_nz = 1, 3, 1
        if not any(f"--{ax}" in " ".join(__import__("sys").argv) for ax in ["nx", "ny", "nz"]):
            args.nx, args.ny, args.nz = default_nx, default_ny, default_nz

        print("=" * 60)
        print("Co2(OH)2CO3: monoclinic P2_1/c → orthorhombic → reoriented")
        print("=" * 60)

        s = build_co2_oh2_co3()
        print_info("1. Monoclinic (original)", s)

        s = monoclinic_to_orthorhombic(s)
        print_info("2. Orthorhombic (beta→90°)", s)

        s = reorient_axes(s)
        print_info("3. Axes reoriented (new_a=old_c, new_b=old_a, new_c=old_b)", s)

        s = make_supercell(s, args.nx, args.ny, args.nz)
        print_info(f"4. Supercell {args.nx}x{args.ny}x{args.nz}", s)

        fname = f"POSCAR_co2_oh2_co3_ortho_{args.nx}x{args.ny}x{args.nz}.vasp"

    s.to(fmt="poscar", filename=fname)
    print(f"  Written to {fname}")


if __name__ == "__main__":
    main()
