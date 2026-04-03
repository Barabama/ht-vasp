"""
HT-VASP - QHA workflow

Quasi-Harmonic Approximation workflows for thermodynamic properties.
"""

import re
import logging
import traceback
from pathlib import Path
from typing import Any, Literal, TypedDict

from pymatgen.core import Structure
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.flows.phonons import PhononMaker
from atomate2.vasp.flows.qha import QhaMaker
from atomate2.vasp.jobs.core import TightRelaxMaker, DielectricMaker
from atomate2.vasp.jobs.phonons import PhononDisplacementMaker
from atomate2.vasp.sets.core import StaticSetGenerator, TightRelaxSetGenerator
from custodian.vasp.handlers import VaspErrorHandler
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
        **kwargs,
    ):
        # Store will be initialized in run_flow to allow custom paths
        self.store = None
        self.worker_name = worker_name

        self.supercell_matrix = supercell_matrix

        # Default INCAR settings
        default_incar = {
            "ENCUT": 450,
            "ISTART": 0,
            "ICHARG": 2,
            # Electronic
            "ISMEAR": 1,
            "SIGMA": 0.1,
            "ALGO": "Normal",
            "NELM": 200,
            "NELMIN": 6,
            "NELMDL": -6,
            # Ionic
            "IBRION": 2,
            "ISIF": 2,
            "NSW": 100,
            "POTIM": 0.2,
            "EDIFF": 1e-6,
            "EDIFFG": -0.01,
            # Magnetic
            "ISPIN": 2,
            # Precision
            "ISYM": 0,
            "LREAL": "Auto",
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
        if global_incar:
            default_incar.update(global_incar)
        global_incar = default_incar
        relax_incar = relax_incar or {}
        eos_incar = eos_incar or {}
        phonon_incar = phonon_incar or {}

        # R3 structural relaxation
        initial_relax_maker = DoubleRelaxMaker.from_relax_maker(
            TightRelaxMaker(
                run_vasp_kwargs={"handlers": [VaspErrorHandler()], **vasp_args},
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=TightRelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
                    user_incar_settings={
                        **global_incar,
                        "NELM": 100,
                        "ISIF": 3,
                        **relax_incar,
                    },
                ),
            )
        )

        # EOS relaxation
        eos_relax_maker = DoubleRelaxMaker.from_relax_maker(
            TightRelaxMaker(
                run_vasp_kwargs={"handlers": [VaspErrorHandler()], **vasp_args},
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=TightRelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
                    user_incar_settings={
                        **global_incar,
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
                    **global_incar,
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
                    **global_incar,
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
            number_of_frames=8,
            ignore_imaginary_modes=True,
        )

    def run_flow(
        self,
        name: str,
        structure: Structure,
        flow_dir: Path | str,
        dir_format: str = "{name}",
        store_path: Path | str = "",
        resume: bool = True,
    ) -> dict[str, Any] | None:
        flow_dir = Path(flow_dir)
        flow_dir.mkdir(parents=True, exist_ok=True)

        # Initialize store
        if not store_path:
            store_path = Path(flow_dir, "store.json").resolve()
        else:
            store_path = Path(store_path).resolve()

        self.store = JobStore(
            JSONStore(str(store_path), read_only=False),
            additional_stores={"data": MemoryStore()},
        )

        flow = self.qha_maker.make(structure, supercell_matrix=self.supercell_matrix)
        log.info(f"Running QHA flow for struct {name} in {flow_dir}")

        try:
            run_locally_custom(
                flow,
                store=self.store,
                root_dir=flow_dir,
                dir_format=dir_format,
                resume=resume,
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

            # Get analyze_free_energy job output
            qha_job_doc = self.store.query_one(
                criteria={"name": {"$regex": "analyze_free_energy"}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},
            )
            if qha_job_doc is None:
                raise ValueError(f"No 'analyze_free_energy' job found in store {store_path}")

            output = self.store.get_output(uuid=qha_job_doc["uuid"], which="last", load=True)
            log.info(f"QHA flow for struct {name} completed successfully")

            output["name"] = name
            output["deformation_energies"] = deformation_energies
            return output

        except Exception as e:
            log.error(f"QHA flow for struct {name} failed: {e}")
            log.error(traceback.format_exc())
            return None

        finally:
            self.close()
