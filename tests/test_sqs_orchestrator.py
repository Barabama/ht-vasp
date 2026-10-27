"""Offline tests for the SQS generation and stage-gating layer."""

from __future__ import annotations

import json
from dataclasses import dataclass

import pytest
from pymatgen.core import Lattice, Structure

from htvasp.model import SqsPlan, SqsSystem, TemplateRecord, TemplateRegistry, TemplateSnapshot, structure_checksum
from htvasp.sqs.generator import GeneratedSqs
from htvasp.sqs.orchestrator import SqsOrchestrationError, SqsOrchestrator, endpoint_cache_key


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


@dataclass
class FakeGenerator:
    engine_name: str = "fake-sqs"

    def generate(
        self,
        template: Structure,
        point,
        *,
        element_a: str,
        element_b: str,
        iterations: int,
        seed: int,
        engine_config=None,
    ) -> GeneratedSqs:
        structure = template.copy()
        for index in range(len(structure)):
            structure.replace(index, element_a if index < point.element_a_count else element_b)
        return GeneratedSqs(
            structure,
            {
                "engine": self.engine_name,
                "engine_version": "test",
                "iterations": iterations,
                "seed": seed,
                "seed_applied": True,
                "objective": 0.0,
            },
        )


def _plans(tmp_path):
    template = _template(tmp_path)
    return [
        SqsPlan.default(SqsSystem.from_name("BCC-Al-Nb"), template, supercell=(2, 2, 2)),
        SqsPlan.default(SqsSystem.from_name("BCC-Al-Ti"), template, supercell=(2, 2, 2)),
    ]


def _observation(spec):
    species = [
        element
        for element, count in spec.structure.composition.items()
        for _ in range(count)
    ]
    final_structure = Structure(
        Lattice.cubic(5.72),
        species,
        [[index / spec.structure.n_sites, 0, 0] for index in range(spec.structure.n_sites)],
    )
    return {
        "calculation_id": spec.calculation_id,
        "state": "successful",
        "energy_total_eV": -100.0,
        "energy_per_atom_eV": -100.0 / spec.structure.n_sites,
        "n_sites": spec.structure.n_sites,
        "final_structure": final_structure.as_dict(),
        "task_document_sha256": "a" * 64,
    }


def test_generate_deduplicates_shared_endpoints_and_writes_manifest(tmp_path):
    orchestrator = SqsOrchestrator(tmp_path / "batch", generator=FakeGenerator())
    manifest = orchestrator.generate(_plans(tmp_path), settings_fingerprint="settings-v1")

    assert manifest.status() == {
        "systems": 2,
        "endpoints": 3,
        "endpoints_successful": 0,
        "endpoints_pending": 3,
        "endpoints_failed": 0,
        "points": 6,
        "points_successful": 0,
        "points_pending": 6,
        "points_failed": 0,
    }
    assert len(orchestrator.pending_endpoints()) == 3
    assert manifest.endpoints[0].consumers or manifest.endpoints[1].consumers
    assert (tmp_path / "batch" / "batch.json").is_file()
    assert len(list((tmp_path / "batch" / "systems").glob("*/manifest.json"))) == 2
    manifest2 = orchestrator.load()
    assert manifest2.fingerprint == manifest.fingerprint


def test_points_are_blocked_until_all_endpoints_succeed(tmp_path):
    orchestrator = SqsOrchestrator(tmp_path / "batch", generator=FakeGenerator())
    manifest = orchestrator.generate(_plans(tmp_path), settings_fingerprint="settings-v1")

    with pytest.raises(SqsOrchestrationError, match="before all endpoints succeed"):
        orchestrator.pending_points()
    with pytest.raises(SqsOrchestrationError, match="does not match requested run"):
        manifest.assert_points_ready(expected_settings_fingerprint="settings-v2")

    for endpoint in manifest.endpoints:
        manifest = orchestrator.mark_endpoint_success(
            endpoint.cache_key,
            _observation(endpoint.spec),
        )

    assert manifest.status()["endpoints_successful"] == 3
    assert len(orchestrator.pending_points()) == 6
    assert all(path.endswith(".observation.json") for path in [
        record.observation_path for record in manifest.endpoints
    ])
    (tmp_path / "batch" / manifest.endpoints[0].observation_path).unlink()
    assert any("missing observation" in error for error in orchestrator.validate())


def test_endpoint_cache_key_changes_with_template_or_settings():
    common = {
        "phase": "BCC",
        "component": "Al",
        "supercell": (2, 2, 2),
    }
    first = endpoint_cache_key(**common, template_checksum="a" * 64, settings_fingerprint="v1")
    second = endpoint_cache_key(**common, template_checksum="b" * 64, settings_fingerprint="v1")
    third = endpoint_cache_key(**common, template_checksum="a" * 64, settings_fingerprint="v2")
    assert len({first, second, third}) == 3


def test_validation_detects_missing_or_corrupt_artifact(tmp_path):
    orchestrator = SqsOrchestrator(tmp_path / "batch", generator=FakeGenerator())
    orchestrator.generate([_plans(tmp_path)[0]], settings_fingerprint="settings-v1")
    endpoint = orchestrator.load().endpoints[0]
    artifact = tmp_path / "batch" / endpoint.spec.structure.path
    artifact.unlink()
    errors = orchestrator.validate()
    assert any("missing structure artifact" in error for error in errors)


def test_generate_requires_explicit_overwrite_for_existing_batch(tmp_path):
    orchestrator = SqsOrchestrator(tmp_path / "batch", generator=FakeGenerator())
    plans = [_plans(tmp_path)[0]]
    orchestrator.generate(plans, settings_fingerprint="settings-v1")
    with pytest.raises(SqsOrchestrationError, match="pass overwrite=True"):
        orchestrator.generate(plans, settings_fingerprint="settings-v1")
    rebuilt = orchestrator.generate(plans, settings_fingerprint="settings-v2", overwrite=True)
    assert rebuilt.endpoints[0].spec.settings_fingerprint == "settings-v2"


def test_collect_recomputes_hmix_and_writes_versioned_dataset(tmp_path):
    orchestrator = SqsOrchestrator(tmp_path / "batch", generator=FakeGenerator())
    manifest = orchestrator.generate(
        [_plans(tmp_path)[0]],
        settings_fingerprint="settings-v1",
        calculation={
            "settings_fingerprint": "settings-v1",
            "workflow": {"worker": "SqsMixingWorker", "flow": "relax_then_static"},
            "relax_incar": {"ISIF": 3},
            "static_incar": {"NSW": 0},
        },
    )
    for endpoint in manifest.endpoints:
        orchestrator.mark_endpoint_success(endpoint.cache_key, _observation(endpoint.spec))
    for point in orchestrator.load().systems[0].points:
        orchestrator.mark_point_success(point.calculation_id, _observation(point.spec))

    output = tmp_path / "mixing.json"
    data = orchestrator.collect("BCC-Al-Nb", output_path=output)
    assert data.schema_version == "htvasp-sqs-mixing/v1"
    assert len(data.points) == 3
    assert data.points[0].h_mix_eV_per_atom == pytest.approx(
        -100.0 / 16 + 100.0 / 16
    )
    assert output.is_file()
    assert type(orchestrator.collect("BCC-Al-Nb", output_path=output)).__name__ == "SqsMixingData"


def test_collect_default_output_matches_endmember_layout(tmp_path):
    root = tmp_path / "BCC-Al-Nb" / "sqsflow"
    orchestrator = SqsOrchestrator(root, generator=FakeGenerator())
    manifest = orchestrator.generate(
        [_plans(tmp_path)[0]],
        settings_fingerprint="settings-v1",
        calculation={
            "settings_fingerprint": "settings-v1",
            "workflow": {"worker": "SqsMixingWorker", "flow": "relax_then_static"},
            "relax_incar": {"ISIF": 3},
            "static_incar": {"NSW": 0},
        },
    )
    for endpoint in manifest.endpoints:
        orchestrator.mark_endpoint_success(endpoint.cache_key, _observation(endpoint.spec))
    for point in orchestrator.load().systems[0].points:
        orchestrator.mark_point_success(point.calculation_id, _observation(point.spec))

    orchestrator.collect("BCC-Al-Nb")
    assert (tmp_path / "BCC-Al-Nb" / "BCC-Al-Nb-sqs.json").is_file()


def test_point_state_and_system_manifest_are_persisted(tmp_path):
    orchestrator = SqsOrchestrator(tmp_path / "batch", generator=FakeGenerator())
    manifest = orchestrator.generate([_plans(tmp_path)[0]], settings_fingerprint="settings-v1")
    for endpoint in manifest.endpoints:
        orchestrator.mark_endpoint_success(endpoint.cache_key, _observation(endpoint.spec))
    point = orchestrator.load().systems[0].points[0]
    orchestrator.mark_point_success(point.calculation_id, _observation(point.spec))

    system_manifest = tmp_path / "batch" / "systems" / "BCC-Al-Nb" / "manifest.json"
    persisted = system_manifest.read_text(encoding="utf-8")
    assert '"state": "successful"' in persisted
    assert orchestrator.load().systems[0].points[0].state == "successful"


def test_failed_point_is_recorded_without_being_collectable(tmp_path):
    orchestrator = SqsOrchestrator(tmp_path / "batch", generator=FakeGenerator())
    manifest = orchestrator.generate([_plans(tmp_path)[0]], settings_fingerprint="settings-v1")
    for endpoint in manifest.endpoints:
        orchestrator.mark_endpoint_success(endpoint.cache_key, _observation(endpoint.spec))
    point = orchestrator.load().systems[0].points[0]
    orchestrator.mark_point_failure(point.calculation_id, "test failure")
    failed = orchestrator.load().systems[0].points[0]
    assert failed.state == "failed"
    assert failed.failure_path is not None
    with pytest.raises(SqsOrchestrationError, match="before all points succeed"):
        orchestrator.collect("BCC-Al-Nb")


def test_collect_rejects_legacy_observations_without_batch_fingerprint(tmp_path):
    orchestrator = SqsOrchestrator(tmp_path / "batch", generator=FakeGenerator())
    manifest = orchestrator.generate([_plans(tmp_path)[0]], settings_fingerprint="settings-v1")
    for endpoint in manifest.endpoints:
        orchestrator.mark_endpoint_success(endpoint.cache_key, _observation(endpoint.spec))
    for point in orchestrator.load().systems[0].points:
        orchestrator.mark_point_success(point.calculation_id, _observation(point.spec))
    persisted = orchestrator.load()
    endpoint_path = tmp_path / "batch" / persisted.endpoints[0].observation_path
    payload = json.loads(endpoint_path.read_text(encoding="utf-8"))
    payload.pop("batch_fingerprint", None)
    endpoint_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(SqsOrchestrationError, match="different or legacy"):
        orchestrator.collect("BCC-Al-Nb")
