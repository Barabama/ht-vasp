"""Atomate2 Maker for one SQS Relax -> Static calculation."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from atomate2.vasp.jobs.core import RelaxMaker, StaticMaker
from atomate2.vasp.sets.core import RelaxSetGenerator, StaticSetGenerator
from jobflow import Flow, Maker
from pymatgen.core import Structure

from htvasp.sqs.jobs import make_static_observation


@dataclass
class SqsCalculationMaker(Maker):
    """Build a ``RelaxMaker -> StaticMaker -> observation`` Flow."""

    name: str = "sqs calculation"
    run_vasp_kwargs: dict[str, Any] = field(default_factory=dict)
    potcar_functional: str = "PBE"
    global_incar: dict[str, Any] = field(default_factory=dict)
    relax_incar: dict[str, Any] = field(default_factory=dict)
    static_incar: dict[str, Any] = field(default_factory=dict)
    relax_reciprocal_density: int | None = None
    static_reciprocal_density: int | None = None

    def __post_init__(self) -> None:
        relax_kpoints = (
            {"reciprocal_density": self.relax_reciprocal_density}
            if self.relax_reciprocal_density is not None
            else {}
        )
        static_kpoints = (
            {"reciprocal_density": self.static_reciprocal_density}
            if self.static_reciprocal_density is not None
            else {}
        )
        self.relax_maker = RelaxMaker(
            name="sqs relax",
            run_vasp_kwargs=self.run_vasp_kwargs,
            stop_children_kwargs={"handle_unsuccessful": False},
            copy_vasp_kwargs={"additional_vasp_files": ("WAVECAR",)},
            input_set_generator=RelaxSetGenerator(
                user_potcar_functional=self.potcar_functional,
                user_kpoints_settings=relax_kpoints,
                user_incar_settings={
                    **self.global_incar,
                    "LWAVE": True,
                    **self.relax_incar,
                    # SQS-mix uses a full cell relaxation in this milestone.
                    "ISIF": 3,
                },
            ),
        )
        self.static_maker = StaticMaker(
            name="sqs static",
            run_vasp_kwargs=self.run_vasp_kwargs,
            stop_children_kwargs={"handle_unsuccessful": False},
            # Static deliberately starts a clean SCF.  Relax changes the
            # lattice, cutoff and k-point mesh, so a Relax WAVECAR is not a
            # reliable restart artifact for this stage.
            copy_vasp_kwargs={},
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=self.potcar_functional,
                user_kpoints_settings=static_kpoints,
                user_incar_settings={
                    **self.global_incar,
                    "ISTART": 0,
                    "ICHARG": 2,
                    "IBRION": -1,
                    "NSW": 0,
                    "ALGO": "Normal",
                    "NELM": 200,
                    "LWAVE": False,
                    "LCHARG": False,
                    "LAECHG": False,
                    **self.static_incar,
                },
            ),
        )

    def make(
        self,
        structure: Structure,
        calculation_spec: Any | None = None,
        *,
        calculation_id: str | None = None,
        prev_dir: str | Path | None = None,
    ) -> Flow:
        """Create the Flow for one endpoint or SQS composition point.

        ``calculation_spec`` is metadata only.  It may be an
        ``SqsCalculationSpec`` or a mapping containing ``calculation_id``.
        For compatibility, a string passed as the second positional argument
        is treated as ``calculation_id``.
        """

        if isinstance(calculation_spec, str) and calculation_id is None:
            calculation_id = calculation_spec
            calculation_spec = None
        spec_id = _calculation_id_from_spec(calculation_spec)
        if calculation_id is not None and spec_id is not None and calculation_id != spec_id:
            raise ValueError(
                "calculation_id does not match calculation_spec.calculation_id"
            )
        resolved_id = calculation_id or spec_id or "sqs-calculation"

        relax_job = self.relax_maker.make(structure, prev_dir=prev_dir)
        static_job = self.static_maker.make(
            relax_job.output.structure,
            prev_dir=relax_job.output.dir_name,
        )
        observation_job = make_static_observation(
            static_job.output,
            calculation_id=resolved_id,
            batch_fingerprint=_spec_metadata(calculation_spec).get("batch_fingerprint"),
            intent_key=_spec_metadata(calculation_spec).get("intent_key"),
        )
        observation_job.name = "sqs observation"
        return Flow(
            [relax_job, static_job, observation_job],
            output=observation_job.output,
            name=self.name,
        )


def _calculation_id_from_spec(spec: Any | None) -> str | None:
    if spec is None:
        return None
    if isinstance(spec, dict):
        value = spec.get("calculation_id")
    else:
        value = getattr(spec, "calculation_id", None)
    if value is None:
        raise ValueError("calculation_spec must provide calculation_id")
    value = str(value).strip()
    if not value:
        raise ValueError("calculation_spec.calculation_id must be non-empty")
    return value


def _spec_metadata(spec: Any | None) -> dict[str, Any]:
    if spec is None:
        return {}
    if isinstance(spec, dict):
        return dict(spec.get("metadata", {}))
    return dict(getattr(spec, "metadata", {}) or {})


__all__ = ["SqsCalculationMaker"]
