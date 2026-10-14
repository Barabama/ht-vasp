"""The endpoint pool: which cells the DFPT-QHA run covers, and why those cells.

Two decisions this module encodes
---------------------------------
**Representation = conventional for everything.**  B2 and L1_2 are simple-cubic
Bravais lattices with a 2- or 4-atom basis, so their cells ARE primitive
already -- `use_symmetrized_structure` is a no-op for all 77 binary endpoints.
The only cells it would inflate are the SER entries whose element is fcc or bcc
(Al, Nb, Ta), and those are inflated on purpose.

**Lattice constant = Vegard from elemental volumes, NOT the `endmember.py`
template.**  The template reuses one lattice constant (AlNi 2.8597 for BCC,
Cu3Au 3.7354 for FCC) for every composition, which is up to -47% in volume.
Vegard is within a few percent.  The occupancy and chemistry still come from
`endmember.py` -- only the lattice constant is overridden, which keeps the
tested generator in the loop while removing the part of it that was wrong.

The duplicate pairs
-------------------
Four binary endpoints are the SAME crystal as an SER entry: an A=A endpoint's
Vegard volume IS the element's own volume, so for Al (fcc) and Cr/Nb/Ta (bcc)
the two cells agree to <0.002% in `a`.  They are one calculation with two
CALPHAD roles, not two independent data points, which is what
`physical_state_id` records.
"""

from __future__ import annotations

import itertools
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

SCHEMA_VERSION = "legacy-qha-pool/pool-structures/v1"

EndpointClass = Literal["BCC_B2", "FCC_L12", "SER"]
PoolPhase = Literal["BCC", "FCC", "SER"]

POOL_ELEMENTS: tuple[str, ...] = ("Al", "Cr", "Hf", "Nb", "Ta", "Ti", "Zr")
ENDPOINT_CLASSES = ("BCC_B2", "FCC_L12", "SER")


@dataclass(frozen=True)
class ElementState:
    """The elemental volume a Vegard estimate is built from."""

    phase: str
    volume_per_atom: float
    mp_id: str


# Ground-state phase and volume per atom, from the MP entries `endmember.py`
# names.  These are the `is_proxy=False` rows.
GROUND_STATE: dict[str, ElementState] = {
    "Al": ElementState("FCC", 16.4717, "mp-134"),
    "Cr": ElementState("BCC", 13.0845, "mp-90"),
    "Hf": ElementState("HCP", 22.2640, "mp-103"),
    "Nb": ElementState("BCC", 18.2581, "mp-75"),
    "Ta": ElementState("BCC", 18.1300, "mp-50"),
    "Ti": ElementState("HCP", 17.3194, "mp-46"),
    "Zr": ElementState("HCP", 23.4997, "mp-131"),
}

# Element volumes in a phase that is NOT its ground state.  Only entries that
# have actually been looked up appear here; anything else falls back to the
# ground-state volume and is flagged `is_proxy` so a reader can see that the
# number is a stand-in rather than a measurement.
POLYMORPH_STATE: dict[tuple[str, str], ElementState] = {
    ("Al", "BCC"): ElementState("BCC", 16.1160, "mp-998860"),
    ("Hf", "BCC"): ElementState("BCC", 22.0870, "mp-100"),
    # fcc polymorphs of the bcc/hcp metals have no looked-up entry yet
}


def volume_per_atom(element: str, phase: str) -> tuple[float, str, bool]:
    """``(volume per atom, provenance, is_proxy)`` for one element in one phase."""

    polymorph = POLYMORPH_STATE.get((element, phase))
    if polymorph is not None:
        return polymorph.volume_per_atom, f"polymorph {polymorph.mp_id}", False
    ground = GROUND_STATE[element]
    if ground.phase == phase:
        return ground.volume_per_atom, f"ground state {ground.mp_id}", False
    return (
        ground.volume_per_atom,
        f"PROXY: ground state {ground.mp_id} is {ground.phase}, used for {phase}",
        True,
    )


def cubic_a(vpa: float, n_atoms: int) -> float:
    """Lattice constant of a cubic cell holding ``n_atoms`` at ``vpa`` each."""

    return (vpa * n_atoms) ** (1.0 / 3.0)


def b2_cell(a: float, element_a: str, element_b: str) -> Any:
    """B2 / A2 conventional cell: element_a at 1a, element_b at 1b."""

    from pymatgen.core import Lattice, Structure

    return Structure(
        Lattice.cubic(a), [element_a, element_b],
        [[0, 0, 0], [0.5, 0.5, 0.5]], coords_are_cartesian=False,
    )


def l12_cell(a: float, element_a: str, element_b: str) -> Any:
    """L1_2 / A1 conventional cell: element_a at 1a, element_b on all of 3c."""

    from pymatgen.core import Lattice, Structure

    return Structure(
        Lattice.cubic(a), [element_a, element_b, element_b, element_b],
        [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
        coords_are_cartesian=False,
    )


def conventional_cell(structure: Any) -> Any:
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    return SpacegroupAnalyzer(structure, symprec=1e-3).get_conventional_standard_structure()


def site_role(structure: Any) -> dict[str, Any]:
    """Which element sits on 1a and which on 3c, read from the fractional coords."""

    import numpy as np

    if len(structure) not in (2, 4):
        return {}
    at_1a = None
    at_3c: list[str] = []
    for site in structure:
        frac = np.array(site.frac_coords) % 1.0
        if np.allclose(frac, 0.0, atol=1e-4) or np.allclose(frac, 1.0, atol=1e-4):
            at_1a = str(site.specie)
        else:
            at_3c.append(str(site.specie))
    return {"site_1a": at_1a, "site_3c": sorted(at_3c) if at_3c else []}


class PoolDecisions(BaseModel):
    """The choices the whole pool rests on, recorded next to the data."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    representation: str
    lattice_constant: str
    template_lattice_used: bool
    duplicate_pairs_share_one_curve: bool


class PoolEndpoint(BaseModel):
    """One endpoint: the cell, how it was built, and what it is."""

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    name: str
    #: serialized as "class" -- that is the key the existing pool files use
    structure_class: EndpointClass = Field(alias="class")
    elements: list[str]
    phase: PoolPhase
    structure: dict[str, Any]
    n_atoms: int = Field(gt=0)
    lattice: str
    a: float
    volume_per_atom: float
    source: str
    vpa_source: dict[str, str] | str | None = None
    is_proxy: bool
    detour: str | None = None
    elem_mpid_note: str | None = None
    site_role: dict[str, Any] | None = None
    #: the crystal this endpoint IS, when it coincides with another endpoint
    physical_state_id: str | None = None

    @property
    def is_duplicate_of_another_endpoint(self) -> bool:
        return self.physical_state_id is not None


class PoolManifest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal[SCHEMA_VERSION] = Field(SCHEMA_VERSION, alias="schema")
    n_endpoints: int
    by_class: dict[str, int]
    decisions: PoolDecisions
    physical_states: dict[str, list[str]]
    duplicate_states: dict[str, list[str]]
    endpoints: dict[str, PoolEndpoint]

    def dump(self) -> dict[str, Any]:
        return self.model_dump(by_alias=True, mode="json")


def enumerate_endpoints() -> list[dict[str, Any]]:
    """The 84 endpoint specifications, before any structure is built."""

    specs: list[dict[str, Any]] = []
    for element_a, element_b in itertools.combinations_with_replacement(POOL_ELEMENTS, 2):
        specs.append({
            "name": f"BCC-{element_a}-{element_b}-2",
            "class": "BCC_B2",
            "elements": [element_a, element_b],
            "phase": "BCC",
        })
    for element_a in POOL_ELEMENTS:
        for element_b in POOL_ELEMENTS:
            specs.append({
                "name": f"FCC-{element_a}-{element_b}-4",
                "class": "FCC_L12",
                "elements": [element_a, element_b],
                "phase": "FCC",
                "site_note": f"{element_a}@1a / {element_b}@3c",
            })
    for element in POOL_ELEMENTS:
        specs.append({
            "name": f"SER-{element}-1",
            "class": "SER",
            "elements": [element],
            "phase": "SER",
        })
    return specs


def physical_state_id(spec: dict[str, Any]) -> str | None:
    """The name of the physical crystal this endpoint IS, for duplicate pairs.

    Only an A=A binary endpoint can coincide with an SER entry, and only when
    the element's ground-state phase IS that phase.  The four that do are Al
    (fcc) and Cr, Nb, Ta (bcc).
    """

    elements = spec["elements"]
    if spec["class"] == "SER":
        return f"{GROUND_STATE[elements[0]].phase.lower()}-{elements[0]}"
    if spec["class"] in ("BCC_B2", "FCC_L12") and elements[0] == elements[1]:
        phase = spec["phase"]
        if GROUND_STATE[elements[0]].phase == phase:
            return f"{phase.lower()}-{elements[0]}"
    return None


def ser_cell(element: str, cache: Path) -> tuple[Any, dict[str, Any]]:
    """The element's MP ground state, fetched through ``endmember.py``.

    The cache is checked first so a rebuild does not need network access.  Any
    element the production ``elem_mpid`` table does not cover is supplied at
    the call site and the fact is recorded -- patching the table from here
    would hide that the default is incomplete.
    """

    from pymatgen.core import Structure

    cached = cache / f"ser_{element}.json"
    if cached.is_file():
        return Structure.from_dict(json.loads(cached.read_text())), {"source": "cache"}

    from htvasp.model.endmember import Endmember

    table = dict(Endmember().elem_mpid)
    note = None
    for key, value in ELEM_MPID_SUPPLEMENT.items():
        if key not in table:
            note = f"supplied {key}={value} at the call site (htvasp default lacks it)"
            table[key] = value
    structure = Endmember(elem_mpid=table).get_ser(element)
    cache.mkdir(parents=True, exist_ok=True)
    cached.write_text(json.dumps(structure.as_dict()))
    return structure, {"source": "endmember.get_ser", "elem_mpid_note": note}


# `endmember.elem_mpid` covers 13 elements and the pool needs 7, of which Hf is
# the one it lacks -- `get_ser("Hf")` raises KeyError.  Supplied here rather
# than patched into the production default, so the gap stays visible.
ELEM_MPID_SUPPLEMENT = {"Hf": "mp-103"}


def build_pool(cache: Path) -> PoolManifest:
    """Build every endpoint cell.  Deterministic given the same ``cache``."""

    from pymatgen.core import Structure

    specs = enumerate_endpoints()
    by_class: dict[str, int] = {}
    for spec in specs:
        by_class[spec["class"]] = by_class.get(spec["class"], 0) + 1
    assert by_class == {"BCC_B2": 28, "FCC_L12": 49, "SER": 7}, by_class

    endpoints: dict[str, PoolEndpoint] = {}
    states: dict[str, list[str]] = {}
    for spec in specs:
        name, endpoint_class = spec["name"], spec["class"]
        record: dict[str, Any]
        if endpoint_class == "SER":
            element = spec["elements"][0]
            raw, source = ser_cell(element, cache)
            cell = conventional_cell(raw)
            record = {
                "structure": cell.as_dict(),
                "n_atoms": len(cell),
                "lattice": "conventional of the MP ground state",
                "a": float(cell.lattice.a),
                "volume_per_atom": float(cell.volume) / len(cell),
                "source": source["source"],
                "elem_mpid_note": source.get("elem_mpid_note"),
                "vpa_source": f"MP ground state {GROUND_STATE[element].mp_id}",
                "is_proxy": False,
                "detour": (
                    f"MP gives {len(raw)} atoms; conventional is {len(cell)}"
                    if len(raw) != len(cell) else None
                ),
            }
        else:
            element_a, element_b = spec["elements"]
            phase = spec["phase"]
            vpa_a, source_a, proxy_a = volume_per_atom(element_a, phase)
            vpa_b, source_b, proxy_b = volume_per_atom(element_b, phase)
            if endpoint_class == "BCC_B2":
                vpa, n_atoms = (vpa_a + vpa_b) / 2.0, 2
                cell = b2_cell(cubic_a(vpa, n_atoms), element_a, element_b)
            else:
                vpa, n_atoms = (vpa_a + 3.0 * vpa_b) / 4.0, 4
                cell = l12_cell(cubic_a(vpa, n_atoms), element_a, element_b)
            record = {
                "structure": cell.as_dict(),
                "n_atoms": n_atoms,
                "lattice": "Vegard (NOT the endmember.py template)",
                "a": float(cell.lattice.a),
                "volume_per_atom": vpa,
                "source": "Vegard from elemental volumes",
                "vpa_source": {element_a: source_a, element_b: source_b},
                "is_proxy": bool(proxy_a or proxy_b),
                "site_role": site_role(cell),
            }
        endpoint = PoolEndpoint.model_validate({
            **record,
            "name": name,
            "class": endpoint_class,
            "elements": list(spec["elements"]),
            "phase": spec["phase"],
            "physical_state_id": physical_state_id(spec),
        })
        endpoints[name] = endpoint
        if endpoint.physical_state_id:
            states.setdefault(endpoint.physical_state_id, []).append(name)

    return PoolManifest(
        n_endpoints=len(endpoints),
        by_class=by_class,
        decisions=PoolDecisions(
            representation="conventional for all 84",
            lattice_constant="Vegard from elemental volumes",
            template_lattice_used=False,
            duplicate_pairs_share_one_curve=True,
        ),
        physical_states=states,
        duplicate_states={k: v for k, v in states.items() if len(v) > 1},
        endpoints=endpoints,
    )


__all__ = [
    "ENDPOINT_CLASSES",
    "ELEM_MPID_SUPPLEMENT",
    "ElementState",
    "EndpointClass",
    "GROUND_STATE",
    "POOL_ELEMENTS",
    "POLYMORPH_STATE",
    "PoolDecisions",
    "PoolEndpoint",
    "PoolManifest",
    "SCHEMA_VERSION",
    "b2_cell",
    "build_pool",
    "conventional_cell",
    "cubic_a",
    "enumerate_endpoints",
    "l12_cell",
    "physical_state_id",
    "ser_cell",
    "site_role",
    "volume_per_atom",
]
