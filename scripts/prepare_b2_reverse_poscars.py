#!/usr/bin/env python3
"""Derive ordered B2 reverse-endmember POSCARs from audited forward POSCARs.

The source POSCAR is read with pymatgen.  Only the chemical species assigned to
fractional (0, 0, 0) and (1/2, 1/2, 1/2) are exchanged; lattice, site order,
coordinates, and site properties are otherwise preserved.  A manifest records
both source and derived content hashes so the generated inputs are auditable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Poscar

ROOT = Path(__file__).resolve().parents[1]
POSCAR_DIR = ROOT / "data" / "poscars"
MANIFEST_PATH = POSCAR_DIR / "b2_reverse_manifest.json"

FORWARD_TO_REVERSE = {
    "BCC-Al-Nb": "BCC-Nb-Al",
    "BCC-Al-Ti": "BCC-Ti-Al",
    "BCC-Al-V": "BCC-V-Al",
    "BCC-Nb-Ti": "BCC-Ti-Nb",
    "BCC-Nb-V": "BCC-V-Nb",
    "BCC-Ti-V": "BCC-V-Ti",
}
ORIGIN = (0.0, 0.0, 0.0)
BODY = (0.5, 0.5, 0.5)
_COORD_TOL = 1e-7


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _coord_key(coords: Any) -> tuple[float, float, float]:
    return tuple(round(float(value) % 1.0, 8) for value in coords)


def _find_site(structure: Structure, target: tuple[float, float, float]) -> int:
    matches = []
    for index, site in enumerate(structure):
        delta = [abs((float(actual) % 1.0) - expected) for actual, expected in zip(site.frac_coords, target)]
        # Periodic coordinates can be represented as 0.999999... rather than 0.
        delta = [min(value, 1.0 - value) for value in delta]
        if max(delta) <= _COORD_TOL:
            matches.append(index)
    if len(matches) != 1:
        raise ValueError(f"expected exactly one site at {target}, found indices {matches}")
    return matches[0]


def _species_name(site: Any) -> str:
    return str(site.specie)


def derive_reverse(source_path: Path, derived_path: Path, reverse_name: str) -> dict[str, Any]:
    if derived_path.exists():
        raise FileExistsError(f"refusing to overwrite existing derived POSCAR: {derived_path}")

    source = Structure.from_file(source_path)
    if len(source) != 2:
        raise ValueError(f"{source_path} must contain exactly two B2 sites, found {len(source)}")
    origin_index = _find_site(source, ORIGIN)
    body_index = _find_site(source, BODY)
    if origin_index == body_index:
        raise ValueError(f"{source_path} maps origin and body-center to the same site")

    source_species = [_species_name(site) for site in source]
    origin_element = source_species[origin_index]
    body_element = source_species[body_index]
    if origin_element == body_element:
        raise ValueError(f"{source_path} is not an heteroatomic B2 endmember")

    reversed_species = list(source_species)
    reversed_species[origin_index], reversed_species[body_index] = (
        reversed_species[body_index],
        reversed_species[origin_index],
    )
    derived = Structure(
        lattice=source.lattice,
        species=reversed_species,
        coords=source.frac_coords,
        coords_are_cartesian=False,
        site_properties=source.site_properties,
    )
    derived_path.parent.mkdir(parents=True, exist_ok=True)
    Poscar(derived, comment=reverse_name, sort_structure=False).write_file(derived_path)

    # Re-read the emitted file to hash and validate the actual VASP artifact.
    emitted = Structure.from_file(derived_path)
    if len(emitted) != len(source) or emitted.lattice != source.lattice:
        raise AssertionError(f"derived structure changed lattice or site count: {derived_path}")
    for old, new in zip(source.frac_coords, emitted.frac_coords):
        if max(abs(float(a) - float(b)) for a, b in zip(old, new)) > _COORD_TOL:
            raise AssertionError(f"derived structure changed fractional coordinates: {derived_path}")
    if [_species_name(site) for site in emitted] != reversed_species:
        raise AssertionError(f"derived POSCAR species order is not the requested site exchange: {derived_path}")

    source_rel = source_path.resolve().relative_to(ROOT.resolve()).as_posix()
    derived_rel = derived_path.resolve().relative_to(ROOT.resolve()).as_posix()
    return {
        "reverse_name": reverse_name,
        "reverse_of": source_path.stem,
        "canonical_phase": "B2_BCC",
        "ordered": True,
        "source_path": source_rel,
        "source_sha256": sha256(source_path),
        "derived_path": derived_rel,
        "derived_sha256": sha256(derived_path),
        "generation_rule": "pymatgen Structure: exchange species at fractional (0,0,0) and (0.5,0.5,0.5); preserve lattice, site order, coordinates, and POTCAR species order",
        "site_mapping": {
            "sublattice_1": {
                "fractional_coordinate": list(ORIGIN),
                "source_element": origin_element,
                "derived_element": body_element,
                "source_site_index": origin_index,
                "derived_site_index": origin_index,
            },
            "sublattice_2": {
                "fractional_coordinate": list(BODY),
                "source_element": body_element,
                "derived_element": origin_element,
                "source_site_index": body_index,
                "derived_site_index": body_index,
            },
        },
        "source_species_order": source_species,
        "derived_species_order": reversed_species,
        "source_fractional_coordinates": [list(map(float, site.frac_coords)) for site in source],
        "derived_fractional_coordinates": [list(map(float, site.frac_coords)) for site in emitted],
        "lattice_matrix": [[float(value) for value in row] for row in source.lattice.matrix],
    }


def build_manifest(*, allow_existing_manifest: bool = False) -> dict[str, Any]:
    if MANIFEST_PATH.exists() and not allow_existing_manifest:
        raise FileExistsError(f"refusing to overwrite existing manifest: {MANIFEST_PATH}")

    records = []
    for forward_name, reverse_name in FORWARD_TO_REVERSE.items():
        source_path = POSCAR_DIR / f"{forward_name}.vasp"
        derived_path = POSCAR_DIR / f"{reverse_name}.vasp"
        if not source_path.is_file():
            raise FileNotFoundError(f"missing forward POSCAR: {source_path}")
        records.append(derive_reverse(source_path, derived_path, reverse_name))

    manifest = {
        "schema_version": "b2-reverse-poscar-manifest-v1",
        "canonical_phase": "B2_BCC",
        "ordered": True,
        "generator": "scripts/prepare_b2_reverse_poscars.py",
        "generation_rule": "pymatgen Structure site-species exchange, not text substitution",
        "targets": records,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-existing-manifest", action="store_true", help="only permit an existing manifest; derived POSCARs are still never overwritten")
    args = parser.parse_args()
    manifest = build_manifest(allow_existing_manifest=args.allow_existing_manifest)
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
