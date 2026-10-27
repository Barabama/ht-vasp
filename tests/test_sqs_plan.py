"""Offline tests for the SQS planning and local-template model layer."""

from __future__ import annotations

import pytest
from pymatgen.core import Lattice, Structure

from htvasp.model import (
    SqsCalculationSpec,
    SqsCompositionPoint,
    SqsPlan,
    SqsStructureArtifact,
    SqsSystem,
    TemplateRecord,
    TemplateRegistry,
    TemplateSnapshot,
    structure_checksum,
)


def _template(tmp_path) -> TemplateSnapshot:
    structure = Structure(
        Lattice.cubic(2.86),
        ["Fe", "Fe"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    path = tmp_path / "BCC.vasp"
    structure.to(filename=path, fmt="poscar")
    record = TemplateRecord(
        phase="BCC",
        path="BCC.vasp",
        source_species=("Fe", "Fe"),
        structure_checksum=structure_checksum(structure),
    )
    return TemplateRegistry(tmp_path, {"BCC": record}).snapshot("BCC")


def _artifact(point: SqsCompositionPoint) -> SqsStructureArtifact:
    return SqsStructureArtifact(
        path=f"structures/{point.label}.vasp",
        structure_checksum="a" * 64,
        n_sites=point.n_sites,
        composition={"Al": point.element_a_count, "Nb": point.element_b_count},
    )


def test_system_names_are_canonicalized_in_both_directions():
    assert SqsSystem.from_name("BCC-Nb-Al").canonical_name == "BCC-Al-Nb"
    assert SqsSystem.from_name("bcc-al-nb").key == ("BCC", "Al", "Nb")


def test_composition_points_require_integer_counts_and_have_stable_labels():
    point = SqsCompositionPoint.from_fraction("1/4", 16)
    assert point.element_a_count == 4
    assert point.element_b_count == 12
    assert point.n_sites == 16
    assert point.label == "x0.25"
    with pytest.raises(ValueError, match="integers"):
        SqsCompositionPoint.from_counts(1.5, 2)
    with pytest.raises(ValueError, match="cannot be represented"):
        SqsCompositionPoint.from_fraction("1/3", 16)


def test_template_registry_validates_local_file_and_structure_checksum(tmp_path):
    snapshot = _template(tmp_path)
    assert snapshot.checksum == snapshot.record.structure_checksum
    assert TemplateRegistry(tmp_path, {"BCC": snapshot.record}).get("BCC") == snapshot.structure

    bad_record = TemplateRecord(
        phase="BCC",
        path="BCC.vasp",
        source_species=("Fe", "Fe"),
        structure_checksum="b" * 64,
    )
    with pytest.raises(ValueError, match="checksum mismatch"):
        TemplateRegistry(tmp_path, {"BCC": bad_record}).snapshot("BCC")
    with pytest.raises(ValueError, match="SHA-256"):
        TemplateRecord("BCC", "BCC.vasp", ("Fe",), structure_checksum="not-a-checksum")


def test_plan_manifest_round_trip_preserves_fingerprint(tmp_path):
    system = SqsSystem.from_name("BCC-Nb-Al")
    plan = SqsPlan.default(system, _template(tmp_path), supercell=(2, 2, 2))
    path = plan.write_manifest(tmp_path / "manifest.json")
    restored = SqsPlan.read_manifest(path)
    assert restored.stable_id == "BCC-Al-Nb"
    assert restored.fingerprint == plan.fingerprint
    assert restored.points == plan.points
    assert restored.template.checksum == plan.template.checksum


def test_endpoint_and_point_specs_have_stable_ids_and_round_trip(tmp_path):
    template = _template(tmp_path)
    system = SqsSystem.from_name("BCC-Al-Nb")
    point = SqsCompositionPoint.from_counts(4, 12)
    endpoint = SqsCalculationSpec.endpoint(
        phase="bcc",
        component="al",
        structure=SqsStructureArtifact(
            path="endpoints/BCC-Al.vasp",
            structure_checksum="b" * 64,
            n_sites=16,
            composition={"Al": 16},
        ),
        settings_fingerprint="settings-v1",
    )
    point_spec = SqsCalculationSpec.point_spec(
        system=system,
        point=point,
        structure=_artifact(point),
        settings_fingerprint="settings-v1",
        supercell=(2, 2, 2),
    )

    assert endpoint.stable_id == "BCC-Al"
    assert point_spec.stable_id == "BCC-Al-Nb-x0.25"
    assert point_spec.fingerprint == SqsCalculationSpec.from_dict(point_spec.as_dict()).fingerprint
    assert SqsCalculationSpec.from_dict(endpoint.as_dict()).stable_id == endpoint.stable_id
    assert template.phase == "BCC"

    with pytest.raises(ValueError, match="calculation_id"):
        SqsCalculationSpec(
            kind="endpoint",
            calculation_id="wrong-id",
            phase="BCC",
            component="Al",
            structure=endpoint.structure,
            settings_fingerprint="settings-v1",
        )
