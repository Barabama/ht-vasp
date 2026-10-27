"""
HT-VASP - QHA workflow

Quasi-Harmonic Approximation workflows for thermodynamic properties.
"""

import copy
import hashlib
import json
import math
import re
import logging
import traceback
from pathlib import Path
from typing import Any, TypedDict

from jobflow import Flow
from pymatgen.core import Structure
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.flows.phonons import PhononMaker
from atomate2.vasp.flows.qha import QhaMaker
from atomate2.vasp.jobs.core import RelaxMaker, DielectricMaker
from atomate2.vasp.jobs.phonons import PhononDisplacementMaker
from atomate2.vasp.sets.core import StaticSetGenerator, RelaxSetGenerator

from htvasp.workflows.base import Worker

log = logging.getLogger(__name__)


class QhaData(TypedDict):
    name: str
    structure: dict[str, Any]
    bulk_modulus: float  # DFT E0 GPa
    volumes: list[float]  # EOS
    temperatures: list[float]  # T K
    thermal_expansion: list[float]  # K^(-1)
    bulk_modulus_temperature: list[float]  # B(T) GPa
    heat_capacity_p_numerical: list[float]  # Cp(T) J/mol/K
    gibbs_temperature: list[float]  # G(T)
    gruneisen_temperature: list[float]  # γ(T)
    volume_temperature: list[float]  # V(T)
    free_energies: list[float]  # F(V,T) kJ/mol
    deformation_energies: list[float]  # E0(V) eV
    total_magnetizations: list[float]  # M(V) µB, sorted by volume ascending
    entropies: list[list[float]]  # S(V,T)
    heat_capacities: list[list[float]]  # Cv(V,T) J/mol/K
    helmholtz_volume: list[list[float]]  # A(V,T)


QHA_SCHEMA_VERSION = "qha-v1"
QHA_PRODUCER_VERSION = "htvasp-qha-v1"
QHA_TEMPERATURE_UNIT = "K"
QHA_GIBBS_UNIT = "eV/reduced_formula_unit"


def _qha_reference_role(structure_name: str) -> str:
    """Map a configured endmember name to an explicit DFT-QHA role."""

    name = structure_name.removesuffix("-qha")
    tokens = name.split("-")
    if len(tokens) == 2 and tokens[0].upper() == "SER" and tokens[1]:
        return "DFT_QHA_SER"
    if len(tokens) == 3 and tokens[0].upper() in {"BCC", "FCC"}:
        phase = tokens[0].upper()
        if tokens[1].upper() == tokens[2].upper():
            return f"DFT_QHA_{phase}_PURE"
        return "DFT_QHA_ORDERED"
    raise ValueError(
        f"cannot assign an explicit QHA reference role from structure name {structure_name!r}"
    )


def _finite_array(value: Any, field: str) -> list[float]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"QHA {field} must be a non-empty list")
    result: list[float] = []
    for index, item in enumerate(value):
        try:
            number = float(item)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"QHA {field}[{index}] must be numeric") from exc
        if not math.isfinite(number):
            raise ValueError(f"QHA {field}[{index}] must be finite")
        result.append(number)
    return result


def canonicalize_qha_result(
    data: dict[str, Any],
    *,
    structure_name: str,
    output_path: Path | str | None = None,
    phase_context: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Attach the production QHA contract before writing the result JSON.

    Current atomate2/phonopy output has one trailing temperature more than the
    derived Gibbs curve.  The production artifact deliberately stores only
    equal-length arrays; the input/output counts and transformation are kept in
    provenance, while the source calculation remains reproducible from the
    job store and workflow metadata.
    """

    if not isinstance(data, dict):
        raise ValueError("QHA result must be a dictionary")
    result = copy.deepcopy(data)
    temperatures = _finite_array(result.get("temperatures"), "temperatures")
    gibbs = _finite_array(result.get("gibbs_temperature"), "gibbs_temperature")
    if any(value < 0.0 for value in temperatures):
        raise ValueError("QHA temperatures must be non-negative Kelvin values")
    if any(right <= left for left, right in zip(temperatures, temperatures[1:])):
        raise ValueError("QHA temperatures must be strictly increasing")
    if len(temperatures) == len(gibbs) + 1:
        input_temperature_count = len(temperatures)
        temperatures = temperatures[: len(gibbs)]
        alignment_operation = "trailing_temperature_dropped"
    elif len(temperatures) == len(gibbs):
        input_temperature_count = len(temperatures)
        alignment_operation = "exact"
    else:
        raise ValueError(
            "QHA temperatures and gibbs_temperature must have equal length "
            "or exactly one additional trailing temperature"
        )
    if len(temperatures) != len(gibbs):  # defensive assertion for future edits
        raise ValueError("QHA production arrays must have equal length")

    role = _qha_reference_role(structure_name)
    effective_phase_context = copy.deepcopy(
        phase_context if phase_context is not None else result.get("phase_context")
    )
    if effective_phase_context is not None and not isinstance(effective_phase_context, dict):
        raise ValueError("QHA phase_context must be a dictionary")
    structure_context = dict(result.get("structure_context") or {})
    structure_context.update(
        {
            "structure_name": structure_name.removesuffix("-qha"),
            "workflow": "qha",
        }
    )
    if effective_phase_context:
        structure_context["phase_context"] = copy.deepcopy(effective_phase_context)
    composition_reduced = result.get("composition_reduced")
    structure = result.get("structure")
    if not isinstance(composition_reduced, dict) and isinstance(structure, dict):
        composition_obj = Structure.from_dict(structure).composition.reduced_composition
        composition_reduced = {
            str(element): float(amount)
            for element, amount in composition_obj.as_dict().items()
        }
    if not isinstance(composition_reduced, dict) or not composition_reduced:
        raise ValueError("QHA production result requires composition_reduced")
    composition_reduced = {
        str(element): float(amount)
        for element, amount in composition_reduced.items()
    }
    formula_atom_count = sum(composition_reduced.values())
    if not math.isfinite(formula_atom_count) or formula_atom_count <= 0:
        raise ValueError("QHA composition_reduced must have a positive atom count")

    result.update(
        {
            "schema_version": QHA_SCHEMA_VERSION,
            "producer": "htvasp.QhaWorker",
            "producer_version": QHA_PRODUCER_VERSION,
            "structure_name": structure_name.removesuffix("-qha"),
            "composition_reduced": composition_reduced,
            "units": {
                "temperature": QHA_TEMPERATURE_UNIT,
                "gibbs_temperature": QHA_GIBBS_UNIT,
            },
            "temperature_alignment": "exact",
            "temperatures": temperatures,
            "gibbs_temperature": gibbs,
            "structure_context": structure_context,
            "reference_context": {
                **dict(result.get("reference_context") or {}),
                "role": role,
                "role_status": "explicit",
                "role_source": "htvasp_workflow_structure_name",
            },
            "@module": result.get("@module", "htvasp.workflows.qha"),
            "@class": result.get("@class", "QhaData"),
            "@version": result.get("@version", QHA_PRODUCER_VERSION),
        }
    )
    if effective_phase_context:
        result["phase_context"] = copy.deepcopy(effective_phase_context)
    if result.get("formula_units") is None and result.get("nsites") is not None:
        result["formula_units"] = float(result["nsites"]) / formula_atom_count

    # This is a digest of the pre-canonical producer payload, not a digest of
    # the final JSON bytes (which would be circular if embedded in the JSON).
    # Inventory/discovery computes the final artifact hash independently.
    producer_payload = copy.deepcopy(data)
    producer_payload_sha256 = hashlib.sha256(
        json.dumps(
            producer_payload,
            sort_keys=True,
            separators=(",", ":"),
            default=str,
        ).encode("utf-8")
    ).hexdigest()
    provenance = dict(result.get("provenance") or {})
    provenance.update(
        {
            "workflow": "qha",
            "schema_version": QHA_SCHEMA_VERSION,
            "source_sha256": producer_payload_sha256,
            "reference_role": role,
            "reference_role_status": "explicit",
            "reference_role_source": "htvasp_workflow_structure_name",
            "temperature_unit": QHA_TEMPERATURE_UNIT,
            "gibbs_temperature_unit": QHA_GIBBS_UNIT,
            "temperature_alignment": "exact",
            "temperature_alignment_operation": alignment_operation,
            "input_temperature_count": input_temperature_count,
            "input_gibbs_count": len(gibbs),
            "output_temperature_count": len(temperatures),
            "output_gibbs_count": len(gibbs),
            "producer_payload_sha256": producer_payload_sha256,
        }
    )
    if effective_phase_context:
        provenance["phase_context"] = copy.deepcopy(effective_phase_context)
    if output_path is not None:
        output = Path(output_path).resolve()
        provenance["output_path"] = str(output)
        provenance["source_path"] = str(output.parent / "qhaflow")
        provenance["source_kind"] = "producer_payload"
    result["provenance"] = provenance
    result["artifact"] = {
        **dict(result.get("artifact") or {}),
        "artifact_id": f"{structure_name.removesuffix('-qha')}:qha",
        "kind": "qha_json",
        "uri": str(Path(output_path).resolve()) if output_path is not None else result.get("structure_name", structure_name),
        "schema_version": QHA_SCHEMA_VERSION,
        "producer": "htvasp.QhaWorker",
        "producer_version": QHA_PRODUCER_VERSION,
        "availability": "present",
    }
    return result


# 原来的E-V.dat, 即DFT静态能量E0(V), 近似于 energy_kJ_mol_t[0]
# 或从store取"phonon static eos deformation *" output["output"]["energy"]


class QhaWorker(Worker):
    """
    Worker for Quasi-Harmonic Approximation calculations.
    """

    def __init__(
        self,
        worker_name: str = "qha-worker",
        vasp_args: dict[str, Any] | None = None,
        potcar_functional="PBE",
        global_incar: dict[str, Any] | None = None,
        relax_incar: dict[str, Any] | None = None,
        eos_incar: dict[str, Any] | None = None,
        phonon_incar: dict[str, Any] | None = None,
        temperature_range: tuple[int, int, int] = (0, 3000, 50),
        supercell_matrix: tuple = ((2, 0, 0), (0, 2, 0), (0, 0, 2)),
        linear_strain: tuple = (-0.05, 0.05),
        number_of_frames: int = 8,
        reference_name: str | None = None,
        phase_context: dict[str, Any] | None = None,
        **kwargs,
    ):
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )

        self.supercell_matrix = supercell_matrix
        self.reference_name = reference_name
        self.phase_context = copy.deepcopy(phase_context) if phase_context is not None else None
        relax_incar = relax_incar or {}
        eos_incar = eos_incar or {}
        phonon_incar = phonon_incar or {}

        # R3 structural relaxation
        initial_relax_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                name="init relax",
                run_vasp_kwargs=self.run_vasp_kwargs,
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=self.potcar_functional,
                    user_incar_settings={
                        **self.global_incar,
                        "ISIF": 3,
                        **relax_incar,
                    },
                ),
            )
        )

        # EOS relaxation
        eos_relax_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                name="eos relax",
                run_vasp_kwargs=self.run_vasp_kwargs,
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=self.potcar_functional,
                    user_incar_settings={
                        **self.global_incar,
                        "ISIF": 2,
                        **eos_incar,
                    },
                ),
            )
        )

        # Phonon displacement maker
        phonon_displacement_maker = PhononDisplacementMaker(
            run_vasp_kwargs=self.run_vasp_kwargs,
            stop_children_kwargs={"handle_unsuccessful": False},
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=self.potcar_functional,
                user_incar_settings={
                    **self.global_incar,
                    "IBRION": -1,
                    "ISIF": 2,
                    "NSW": 0,
                    "EDIFF": 1e-7,
                    **phonon_incar,
                },
            ),
        )

        # Phonon maker
        phonon_maker = PhononMaker(
            bulk_relax_maker=None,
            born_maker=None,
            static_energy_maker=phonon_displacement_maker,
            phonon_displacement_maker=phonon_displacement_maker,
            use_symmetrized_structure="conventional",
            sym_reduce=False,
            generate_frequencies_eigenvectors_kwargs={
                "tmin": temperature_range[0],
                "tmax": temperature_range[1],
                "tstep": temperature_range[2],
            },
        )

        # Dielectric maker
        dielectric_maker = DielectricMaker(
            run_vasp_kwargs=self.run_vasp_kwargs,
            stop_children_kwargs={"handle_unsuccessful": False},
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=self.potcar_functional,
                user_incar_settings={
                    **self.global_incar,
                    "IBRION": 6,
                    "NELM": 200,
                    "ISIF": 3,
                    "NSW": 0,
                    "NFREE": 2,
                    "LWAVE": False,
                    "LCHARG": False,
                },
                lepsilon=True,
            ),
        )

        self.flow_maker = QhaMaker(
            initial_relax_maker=initial_relax_maker,
            eos_relax_maker=eos_relax_maker,
            phonon_maker=phonon_maker,
            linear_strain=linear_strain,
            number_of_frames=number_of_frames,
            ignore_imaginary_modes=True,
            min_length=None,
        )

    def _make_flow(self, structure: Structure, prev_dir: Path | str | None = None) -> Flow:
        return self.flow_maker.make(
            structure,
            supercell_matrix=self.supercell_matrix,
            prev_dir=prev_dir,
        )

    def get_result(self, output_job_name: str = "analyze_free_energy") -> dict[str, Any] | None:
        """Override to include deformation energies & magnetizations, sorted by volume ascending."""
        result = self._query_store(output_job_name)
        if result is not None:
            en, mag = self._get_deformation_data()
            result["deformation_energies"] = en
            result["total_magnetizations"] = list(mag) if mag else []
        self.close()
        return result

    def write_result(self, data: dict | None, json_path: Path | str):
        """Write a canonical, production-contract QHA JSON artifact."""

        if not data:
            log.warning("No data to write")
            return
        destination = Path(json_path).resolve()
        structure_name = self.reference_name or destination.stem.removesuffix("-qha")
        canonical = canonicalize_qha_result(
            data,
            structure_name=structure_name,
            output_path=destination,
            phase_context=self.phase_context,
        )
        super().write_result(data=canonical, json_path=destination)

    def _get_deformation_data(self) -> tuple[list[float], list[float]]:
        """
        Extract deformation energies and total magnetizations from the store,
        sorted by volume ascending.

        Returns (energies_eV, magnetizations_uB) matching the order of
        ``analyze_free_energy``'s ``volumes``.
        """
        try:
            if not self.store:
                raise ValueError("Store is not initialized. Run the flow first.")

            entries: list[tuple[int, float, float, float]] = []  # (deform_index, volume, energy, mag)
            for doc in self.store.query(
                criteria={"name": {"$regex": r"phonon static eos deformation \d+"}},
                properties=["uuid", "name"],
            ):
                uuid = doc.get("uuid", "N/A")
                name = doc.get("name", "")
                match = re.search(r"deformation (\d+)", name)
                if not match:
                    log.warning(f"Could not extract deformation index from job name: {name}")
                    continue
                output = self.store.get_output(uuid=uuid, which="last", load=True)
                volume = output.get("volume", 0.0)
                energy = output["output"]["energy"]

                calcs = output.get("calcs_reversed", [])
                if calcs:
                    oucar = calcs[0].get("output", {}).get("outcar", {})
                    if isinstance(oucar, dict):
                        mag = oucar.get("total_magnetization", None)
                    else:
                        mag = None
                else:
                    mag = None

                entries.append((int(match.group(1)), volume, energy, mag))

            # Sort by volume ascending, matching the order used by analyze_free_energy
            entries.sort(key=lambda x: x[1])
            energies = [e for _, _, e, _ in entries]
            magnetizations = [m for _, _, _, m in entries]
            return energies, magnetizations
        except Exception as e:
            log.error(f"Failed to get deformation data: {e}")
            log.error(traceback.format_exc())
            return [], []
