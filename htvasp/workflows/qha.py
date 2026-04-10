"""
HT-VASP - QHA workflow

Quasi-Harmonic Approximation workflows for thermodynamic properties.
"""

import re
import logging
import traceback
from pathlib import Path
from typing import Any, Literal, TypedDict

from jobflow import Flow
from pymatgen.core import Structure
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.flows.phonons import PhononMaker
from atomate2.vasp.flows.qha import QhaMaker
from atomate2.vasp.jobs.core import RelaxMaker, DielectricMaker
from atomate2.vasp.jobs.phonons import PhononDisplacementMaker
from atomate2.vasp.sets.core import StaticSetGenerator, RelaxSetGenerator
from custodian.vasp.handlers import VaspErrorHandler

from htvasp.workflows.base import Worker

log = logging.getLogger(__name__)


class QhaData(TypedDict):
    name: str
    structure: dict[str, Any]
    bulk_modulus: float  # DFT E0 GPa
    volumes: float  # EOS
    temperatures: list[float]  # T K
    thermal_expansion: list[float]  # K^(-1)
    bulk_modulus_temperature: list[float]  # B(T) GPa
    heat_capacity_p_numerical: list[float]  # Cp(T) J/mol/K
    gibbs_temperature: list[float]  # G(T)
    gruneisen_temperature: list[float]  # γ(T)
    volume_temperature: list[float]  # V(T)
    free_energies: list[float]  # F(V,T) kJ/mol
    deformation_energies: list[float]  # E0(T) kJ/mol
    entropies: list[list[float]]  # S(V,T)
    heat_capacities: list[list[float]]  # Cv(V,T) J/mol/K
    helmholtz_volume: list[list[float]]  # A(V,T)


# 原来的E-V.dat, 即DFT静态能量E0(V), 近似于 energy_kJ_mol_t[0]
# 或从store取"phonon static eos deformation *" output["output"]["energy"]


class QhaWorker(Worker):
    """
    Worker for Quasi-Harmonic Approximation calculations.
    """

    def __init__(
        self,
        worker_name: str,
        vasp_args: dict[str, Any],
        potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64",
        global_incar: dict[str, Any] | None = None,
        relax_incar: dict[str, Any] | None = None,
        eos_incar: dict[str, Any] | None = None,
        phonon_incar: dict[str, Any] | None = None,
        temperature_range: tuple[int, int, int] = (0, 3000, 50),
        supercell_matrix: tuple = ((2, 0, 0), (0, 2, 0), (0, 0, 2)),
        execution_mode: Literal["local", "fireworks"] = "local",
        fireworks_config_dir: Path | str | None = None,
        **kwargs,
    ):
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
            execution_mode=execution_mode,
            fireworks_config_dir=fireworks_config_dir,
        )

        self.supercell_matrix = supercell_matrix
        relax_incar = relax_incar or {}
        eos_incar = eos_incar or {}
        phonon_incar = phonon_incar or {}

        # R3 structural relaxation
        initial_relax_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                run_vasp_kwargs={"handlers": [VaspErrorHandler()], **vasp_args},
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
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
                run_vasp_kwargs={"handlers": [VaspErrorHandler()], **vasp_args},
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
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
            run_vasp_kwargs={"handlers": [VaspErrorHandler()], **vasp_args},
            stop_children_kwargs={"handle_unsuccessful": False},
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=potcar_functional,
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
            run_vasp_kwargs={"handlers": [VaspErrorHandler()], **vasp_args},
            stop_children_kwargs={"handle_unsuccessful": False},
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=potcar_functional,
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
            min_length=None,
            number_of_frames=8,
            ignore_imaginary_modes=True,
        )

    def _make_flow(self, structure: Structure) -> Flow:
        return self.flow_maker.make(structure, supercell_matrix=self.supercell_matrix)

    def get_result(self, output_job_name: str = "analyze_free_energy") -> dict[str, Any] | None:
        try:
            output = super().get_result(output_job_name)
            if not output:
                raise ValueError(f"Flow return None for job {output_job_name}")
            if not self.store:
                raise ValueError("Store is not initialized. Run the flow first.")
            self.store.connect()

            # Get all "phonon static eos deformation *" jobs and sort by deformation index
            deformation_data = []
            for doc in self.store.query(
                criteria={"name": {"$regex": r"phonon static eos deformation \d+"}},
                properties=["uuid", "index", "name"],
            ):
                uuid = doc.get("uuid", "N/A")
                name = doc.get("name", "")
                match = re.search(r"deformation (\d+)", name)
                if not match:
                    log.warning(f"Could not extract deformation index from job name: {name}")
                    continue
                output = self.store.get_output(uuid=uuid, which="last", load=True)
                deformation_data.append((int(match.group(1)), output["output"]["energy"]))
            deformation_data.sort(key=lambda x: x[0])
            deformation_energies = [item[1] for item in deformation_data]

            output["deformation_energies"] = deformation_energies
            return output

        except Exception as e:
            log.error(f"Failed to get result for job: {e}")
            log.error(traceback.format_exc())
            return None
        finally:
            self.close()
