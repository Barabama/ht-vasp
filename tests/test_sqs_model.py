"""Offline tests for the SQS model and htvasp-sqs-mixing/v1 contract."""

from __future__ import annotations

from htvasp.model import CalculationProvenance, SqsMixingData, SqsPlan, SqsSystem


def _dataset() -> SqsMixingData:
    plan = SqsPlan.default(SqsSystem.from_name("BCC-Al-Nb"))
    references = {"Al": {"energy": -16.0}, "Nb": {"energy": -32.0}}
    points = {
        point.label: {
            "energy": -point.element_a_count - 2 * point.element_b_count + 0.8 * point.n_sites
        }
        for point in plan.points
    }
    calculation = CalculationProvenance(
        settings_fingerprint="settings-123",
        workflow={"flow": "single_relax_then_static"},
        relax_incar={},
        static_incar={},
    )
    structures = {point.label: plan.expanded_template() for point in plan.points}
    return SqsMixingData.from_runs(
        plan,
        references,
        points,
        point_structures=structures,
        calculation=calculation,
    )


def test_sqs_data_round_trip_has_explicit_calphad_mapping(tmp_path):
    data = _dataset()
    assert data.schema_version == "htvasp-sqs-mixing/v1"
    assert data.calphad_mapping.phase_name == "BCC_A2"
    assert len(data.points) == 3
    assert [point.x_A for point in data.points] == [0.25, 0.5, 0.75]
    path = data.write_json(tmp_path / "al-nb.json")
    restored = SqsMixingData.read_json(path)
    assert restored.calphad_mapping.components == ("AL", "NB")
    assert restored.points[2].static.state == "successful"
    assert len(restored.points[0].sqs.structure_sha256) == 64


def test_sqs_names_are_canonicalized():
    assert SqsSystem.from_name("BCC-Nb-Al").canonical_name == "BCC-Al-Nb"
    assert SqsSystem.from_name("BCC-Nb-Al").key == ("BCC", "Al", "Nb")


def test_default_supercells_reflect_template_site_counts():
    expected = {
        "BCC": ((2, 2, 2), 16),
        "FCC": ((2, 2, 2), 32),
        "HCP": ((2, 2, 1), 32),
    }
    for phase, (supercell, n_sites) in expected.items():
        plan = SqsPlan.default(SqsSystem.from_name(f"{phase}-Al-Nb"))
        assert plan.supercell == supercell
        assert plan.n_sites == n_sites
