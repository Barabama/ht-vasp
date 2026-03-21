"""
HT-VASP - QHA workflow

Quasi-Harmonic Approximation workflows for thermodynamic properties.
"""

import logging
import traceback
from pathlib import Path
from typing import Any, Literal, TypedDict

from pymatgen.core import Structure
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.flows.phonons import PhononMaker
from atomate2.vasp.flows.qha import QhaMaker
from atomate2.vasp.jobs.core import RelaxMaker, DielectricMaker
from atomate2.vasp.jobs.phonons import PhononDisplacementMaker
from atomate2.vasp.sets.core import StaticSetGenerator, RelaxSetGenerator
from maggma.stores import JSONStore, MemoryStore
from jobflow.core.store import JobStore

from htvasp.workflows.base import Worker
from htvasp.utils.run_locally import run_locally_custom

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


class QHAWorker(Worker):
    """
    Worker for Quasi-Harmonic Approximation calculations.
    """

    def __init__(
        self,
        worker_name: str,
        vasp_args: dict[str, Any],
        potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_54",
        incar_settings: dict[str, Any] | None = None,
        temperature_range: tuple[int, int, int] = (0, 3000, 100),
        supercell_matrix: tuple = ((2, 0, 0), (0, 2, 0), (0, 0, 2)),
        **kwargs,
    ):
        # Store will be initialized in run_relax to allow custom paths
        self.store = None

        self.supercell_matrix = supercell_matrix

        # Default INCAR settings
        default_incar = {
            "ENCUT": 400,
            "ISTART": 0,
            "ICHARG": 2,
            # Electronic
            "ISMEAR": 1,
            "SIGMA": 0.2,
            "ALGO": "Fast",
            "NELM": 100,
            "NELMIN": 4,
            "NELMDL": -12,
            # Ionic
            "IBRION": 2,
            "ISIF": 3,
            "NSW": 10,
            "POTIM": 0.2,
            "EDIFF": 1e-5,
            "EDIFFG": -0.02,
            # Magnetic
            "ISPIN": 2,
            "AMIX": 0.04,
            "BMIX": 1e-4,
            "AMIX_MAG": 0.8,
            "BMIX_MAG": 1e-4,
            # Precision
            "KPAR": 2,
            "NCORE": 1,
            "ISYM": 2,
            "LREAL": False,
            "PREC": "Accurate",
            "SYMPREC": 1e-5,
            # Output
            "LWAVE": False,
            "LCHARG": False,
            "LORBIT": None,
            "LOPTICS": False,
            "LVTOT": False,
        }
        # Override with custom settings
        if incar_settings:
            default_incar.update(incar_settings)
        incar_settings = default_incar

        # R3 structural relaxation
        initial_relax_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                run_vasp_kwargs=vasp_args,
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
                    user_incar_settings={
                        **incar_settings,
                        "ALGO": "Fast",
                        "IBRION": 2,
                        "NELM": 200,
                        "ISIF": 3,
                        "NSW": 20,
                        "EDIFF": 1e-7,
                        "EDIFFG": -0.01,
                    },
                ),
            )
        )

        # EOS relaxation
        eos_relax_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                run_vasp_kwargs=vasp_args,
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
                    user_incar_settings={
                        **incar_settings,
                        "ALGO": "Normal",
                        "IBRION": 2,
                        "NELM": 200,
                        "ISIF": 2,
                        "NSW": 0,
                        "EDIFF": 1e-5,
                        "EDIFFG": -0.01,
                    },
                ),
            )
        )

        # Phonon displacement maker
        phonon_displacement_maker = PhononDisplacementMaker(
            run_vasp_kwargs=vasp_args,
            stop_children_kwargs={"handle_unsuccessful": False},
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=potcar_functional,
                user_incar_settings={
                    **incar_settings,
                    "ALGO": "Normal",
                    "IBRION": -1,
                    "NELM": 200,
                    "ISIF": 3,
                    "NSW": 0,
                    "EDIFF": 1e-7,
                    "EDIFFG": -0.01,
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
            run_vasp_kwargs=vasp_args,
            stop_children_kwargs={"handle_unsuccessful": False},
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=potcar_functional,
                user_incar_settings={
                    **incar_settings,
                    "ALGO": "Normal",
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

        # QHA maker
        self.qha_maker = QhaMaker(
            initial_relax_maker=initial_relax_maker,
            eos_relax_maker=eos_relax_maker,
            phonon_maker=phonon_maker,
            min_length=None,
            ignore_imaginary_modes=True,
        )

    def run_flow(
        self,
        struct_name: str,
        structure: Structure,
        flow_dir: Path | str,
        dir_format: str = "{name}",
        store_path: Path | str = "",
    ) -> dict[str, Any] | None:
        flow_dir = Path(flow_dir)
        flow_dir.mkdir(parents=True, exist_ok=True)

        # Initialize store
        if not store_path:
            store_path = Path(flow_dir, "store.json").resolve()
        else:
            store_path = Path(store_path).resolve()
        self.store = JobStore(
            JSONStore(store_path, read_only=False),
            additional_stores={"data": MemoryStore()},
        )

        flow = self.qha_maker.make(structure, supercell_matrix=self.supercell_matrix)
        log.info(f"Running QHA flow for structure {struct_name} in {flow_dir}")

        try:
            run_locally_custom(
                flow,
                store=self.store,
                root_dir=flow_dir,
                dir_format=dir_format,
            )

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

            job = self.store.query_one(
                criteria={"name": {"$regex": "analyze_free_energy"}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},
            )
            if job is None:
                raise ValueError(f"No 'analyze_free_energy' job found in store {store_path}")

            output = self.store.get_output(uuid=job["uuid"], which="last", load=True)
            log.info(f"QHA flow for structure {struct_name} completed successfully")

            return QhaData(
                name=struct_name,
                structure=output["structure"],
                bulk_modulus=output["bulk_modulus"],
                volumes=output["volumes"],
                temperatures=output["temperatures"],
                thermal_expansion=output["thermal_expansion"],
                bulk_modulus_temperature=output["bulk_modulus_temperature"],
                heat_capacity_p_numerical=output["heat_capacity_p_numerical"],
                gibbs_temperature=output["gibbs_temperature"],
                gruneisen_temperature=output["gruneisen_temperature"],
                volume_temperature=output["volume_temperature"],
                free_energies=output["free_energies"],
                deformation_energies=deformation_energies,
                entropies=output["entropies"],
                heat_capacities=output["heat_capacities"],
                helmholtz_volume=output["helmholtz_volume"],
            )

        except Exception as e:
            log.error(f"QHA flow for structure {struct_name} failed: {e}")
            log.error(traceback.format_exc())
            return None

        finally:
            self.close()
