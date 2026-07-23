"""
HT-VASP - QHA workflow

Quasi-Harmonic Approximation workflows for thermodynamic properties.
"""

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
        potcar_functional="PBE_64",
        global_incar: dict[str, Any] | None = None,
        relax_incar: dict[str, Any] | None = None,
        eos_incar: dict[str, Any] | None = None,
        phonon_incar: dict[str, Any] | None = None,
        temperature_range: tuple[int, int, int] = (0, 3000, 50),
        supercell_matrix: tuple = ((2, 0, 0), (0, 2, 0), (0, 0, 2)),
        linear_strain: tuple = (-0.05, 0.05),
        number_of_frames: int = 8,
        **kwargs,
    ):
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )

        self.supercell_matrix = supercell_matrix
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
