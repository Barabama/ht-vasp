# %%
# src/workflow/flow_qha.py

import re
import logging
import traceback
from pathlib import Path
from typing import Any, TypedDict

from pymatgen.core import Structure
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.flows.phonons import PhononMaker
from atomate2.vasp.flows.qha import QhaMaker
from atomate2.vasp.jobs.core import RelaxMaker, DielectricMaker
from atomate2.vasp.jobs.phonons import PhononDisplacementMaker
from atomate2.vasp.sets.core import StaticSetGenerator, RelaxSetGenerator
from maggma.stores import MongoStore, MemoryStore
from jobflow.core.store import JobStore
from jobflow.managers.local import run_locally

log = logging.getLogger(__name__)


class QhaData(TypedDict):
    name: str
    structure: dict[str, Any]
    bulk_modulus_dft_GPa: float  # DFT E0 GPa
    volume_EOS: float  # EOS
    temperatures_K: list[float]  # T K
    thermal_expansion_t: list[float]  # K^(-1)
    bulk_modulus_GPa_t: list[float]  # B(T) GPa
    heat_capp_J_mol_K_t: list[float]  # Cp(T) J/mol/K
    gibbs_t: list[float]  # G(T)
    gruneisen_t: list[float]  # γ(T)
    volume_t: list[float]  # V(T)
    energy_kJ_mol_t: list[float]  # F(V,T) kJ/mol
    energy0_kJ_mol_t: list[float]  # E0(T) kJ/mol
    entropy_v_t: list[list[float]]  # S(V,T)
    heat_capp_J_mol_K_v_t: list[list[float]]  # Cv(V,T) J/mol/K
    helmholtz_v_t: list[list[float]]  # A(V,T)


# 原来的E-V.dat, 即DFT静态能量E0(V), 近似于 energy_kJ_mol_t[0]
# 或从store取"phonon static eos deformation *" output["output"]["energy"]


class QhaWorker:
    def __init__(
        self,
        vasp_args: dict[str, Any],
        potcar_functional: Any,
        incar_settings: dict[str, Any] = {},
        temperature_range: tuple[int, int, int] = (0, 3000, 100),
        supercell_matrix: tuple = ((2, 0, 0), (0, 2, 0), (0, 0, 2)),
    ):
        self.supercell_matrix = supercell_matrix
        incar_settings = {
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
            # Overrides
            **incar_settings,
        }
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
        # Phonon maker
        phonon_maker = PhononMaker(
            bulk_relax_maker=None,
            born_maker=None,  # dielectric_maker,
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
        # QHA maker
        self.qha_maker = QhaMaker(
            initial_relax_maker=initial_relax_maker,
            eos_relax_maker=eos_relax_maker,
            phonon_maker=phonon_maker,
            # linear_strain=(-0.05, 0.05),
            # number_of_frames=4,
            min_length=None,
            ignore_imaginary_modes=True,
        )
        self.store = JobStore(MemoryStore(), additional_stores={"data": MemoryStore()})

    def run_qha(
        self,
        sys_name: str,
        struct: Structure,
        flowdir: Path,
    ) -> QhaData | None:

        flow = self.qha_maker.make(struct, supercell_matrix=self.supercell_matrix)
        log.info(f"System {sys_name} running qha flow")
        try:
            run_locally(
                flow,
                store=self.store,
                create_folders=True,
                root_dir=flowdir,
                # ensure_success=True,
                # raise_immediately=True,
            )

            self.store.connect()

            # Debug: list all jobs (without 'state')
            log.info("\nAll jobs in store:")
            for doc in self.store.query(properties=["name", "uuid", "index"]):
                name = doc.get("name", "N/A")
                uuid = doc.get("uuid", "N/A")
                index = doc.get("index", "N/A")
                log.info(f"  Name: {name} | Index: {index} | UUID: {uuid}")

            
            # Get all "phonon static eos deformation *" jobs and sort by deformation index
            deformation_data = []
            for doc in self.store.query(
                criteria={"name": {"$regex": r"phonon static eos deformation \d+"}},
                properties=["uuid", "name"]
            ):
                uuid = doc.get("uuid", "N/A")
                name = doc.get("name", "")
                # e.g., "phonon static eos deformation 1"
                match = re.search(r'deformation (\d+)', name)
                if not match:
                    log.warning(f"Could not extract deformation index from job name: {name}")
                    continue
                output = self.store.get_output(uuid=uuid, which="last", load=True)
                energy = output["output"]["energy"]
                idx = int(match.group(1))
                deformation_data.append((idx, energy))

            deformation_data.sort(key=lambda x: x[0])
            energy0_list = [item[1] for item in deformation_data]

            # energy0_list = []
            # inputs = ""
            # while inputs != "q":
            #     inputs = input("Enter UUID or q to quit: ")
            #     if inputs == "q":
            #         break
            #     log.info(f"System {sys_name} querying job with UUID: {inputs}")
            #     job = self.store.query_one(
            #         criteria={"uuid": inputs.strip()},
            #         properties=["uuid", "index", "name"],
            #         sort={"index": -1},  # latest index first
            #     )
            #     if job is None:
            #         log.error(f"System {sys_name} no job found with UUID: {inputs}")
            #         continue
            #     output = self.store.get_output(uuid=job["uuid"], which="last", load=True)
            #     with open("job_output.json", "w", encoding="utf-8") as jf:
            #         json.dump(output, jf, indent=2, cls=DateTimeEncoder)
            #     # print(f"System {sys_name} job output: {output}")
            #     # log.info(f"System {sys_name} job output: {output}")
            #     if "output" in output and "energy" in output["output"]:
            #         energy0_list.append(output["output"]["energy"])

            job = self.store.query_one(
                criteria={"name": {"$regex": "analyze_free_energy"}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},  # latest index first
            )
            if job is None:
                raise ValueError(f"System {sys_name} no 'analyze_free_energy' job found")
            output = self.store.get_output(uuid=job["uuid"], which="last", load=True)
            log.info(f"System {sys_name} Qha flow running successfully")
            self.store.close()

            return QhaData(
                name=sys_name,
                structure=output["structure"],
                bulk_modulus_dft_GPa=output["bulk_modulus"],
                volume_EOS=output["volumes"],
                temperatures_K=output["temperatures"],
                thermal_expansion_t=output["thermal_expansion"],
                bulk_modulus_GPa_t=output["bulk_modulus_temperature"],
                heat_capp_J_mol_K_t=output["heat_capacity_p_numerical"],
                gibbs_t=output["gibbs_temperature"],
                gruneisen_t=output["gruneisen_temperature"],
                volume_t=output["volume_temperature"],
                energy_kJ_mol_t=output["free_energies"],
                energy0_kJ_mol_t=energy0_list,
                entropy_v_t=output["entropies"],
                heat_capp_J_mol_K_v_t=output["heat_capacities"],
                helmholtz_v_t=output["helmholtz_volume"],
            )
        except Exception as e:
            log.error(f"System {sys_name} error running qha flow: {e}")
            log.error(traceback.format_exc())
            return None


# %%
# main_qha.py

import json
import shutil
import logging
from pathlib import Path
from datetime import datetime

from custodian.vasp.handlers import VaspErrorHandler
from mp_api.client import MPRester
from pymatgen.core import Element, Structure

# from src.workflow.flow_qha import QhaData, QhaWorker

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


class Modeler:
    def __init__(self, api_key: str):
        self.api_key = api_key
        self.elem_mpid = {
            "Co": "mp-54",  # HCP
            "Cr": "mp-90",  # BCC
            "Fe": "mp-13",  # BCC
            "Mn": "mp-1055908",  # BCC
            "Ni": "mp-23",  # FCC
        }
        self.phase_func = {
            "SER": self.get_ser,
            "BCC": self.get_bcc,
            "FCC": self.get_fcc,
            "HCP": self.get_hcp,
        }

    def get_structure(self, mpid: str) -> Structure:
        """Fetches the structure from the Materials Project with mpid.

        Args:
            mpid: Materials Project ID
        Returns:
            Structure
        """
        with MPRester(self.api_key) as mpr:
            structure = mpr.get_structure_by_material_id(mpid)
        return structure

    def get_ser(self, element: str) -> Structure:
        """
        Fetches the structure SER of element and mpid.

        Args:
            element: Element to fetch structure for
        Returns:
            Structure SER
        """
        mpid = self.elem_mpid[element]
        struct = self.get_structure(mpid)
        log.info(f"Fetched SER for {element}: {mpid}")
        return struct

    def get_bcc(self, elem1: str, elem2: str) -> Structure:
        """
        Creates a BCC alloy structure from two elements.
        Template is Ni_1a-Al_1b, Pm-3m, mp-1487.

        Args:
            elem1: First element symbol
            elem2: Second element symbol
        Returns:
            Structure
        """
        templ = self.get_structure("mp-1487")
        struct = templ.replace_species({"Ni": Element(elem1), "Al": Element(elem2)})
        return struct.get_sorted_structure(lambda site: str(site.specie).lower())

    def get_fcc(self, elem1: str, elem2: str) -> Structure:
        """
        Creates a FCC alloy structure from two elements.
        Template is Au_1a-Cu_3c, Pm-3m, mp-2258.

        Args:
            elem1: First element symbol
            elem2: Second element symbol
        Returns:
            Structure
        """
        templ = self.get_structure("mp-2258")
        struct = templ.replace_species({"Au": Element(elem1), "Cu": Element(elem2)})
        return struct.get_sorted_structure(lambda site: str(site.specie).lower())

    def get_hcp(self, elem1: str, elem2: str) -> Structure:
        """
        Creates a HCP alloy structure from two elements.
        Template is Sn_2c-Ni_6h, P63/mmc, mp-20112.

        Args:
            elem1: First element symbol
            elem2: Second element symbol
        Returns:
            Structure
        """
        templ = self.get_structure("mp-20112")
        struct = templ.replace_species({"Sn": Element(elem1), "Ni": Element(elem2)})
        return struct.get_sorted_structure(lambda site: str(site.specie).lower())

    def get_poscar(self, name: str, outdir: Path | str) -> Structure:
        """
        Reads a POSCAR file and returns the structure.

        Args:
            name: Name of the POSCAR file
            outdir: Directory where the POSCAR file is located
        Returns:
            Structure
        """
        outdir = Path(outdir) if isinstance(outdir, str) else outdir
        poscar = outdir.joinpath(f"{name}.vasp")
        if poscar.is_file():
            return Structure.from_file(poscar)
        else:
            parts = name.split("-")
            phase = parts[0]
            elems = parts[1:]
            struct = self.phase_func[phase](*elems)
            struct.to(fmt="poscar", filename=poscar)
            return struct


def main():
    api_key = "bqqHJQWs8wPyDZnrTgtsWNLevLmIq4MU"

    # names = [
    #     # "BCC-Fe-Fe",
    #     "BCC-Fe-Mn",
    #     "BCC-Fe-Ni",
    #     "BCC-Mn-Mn",
    #     "BCC-Mn-Ni",
    #     # "BCC-Ni-Ni",
    #     # "FCC-Fe-Fe",
    #     "FCC-Fe-Mn",
    #     "FCC-Fe-Ni",
    #     "FCC-Mn-Fe",
    #     "FCC-Mn-Mn",
    #     "FCC-Mn-Ni",
    #     "FCC-Ni-Fe",
    #     "FCC-Ni-Mn",
    #     # "FCC-Ni-Ni",
    #     # "HCP-Fe-Fe",
    #     # "HCP-Fe-Mn",
    #     # "HCP-Fe-Ni",
    #     # "HCP-Mn-Fe",
    #     "HCP-Mn-Mn",
    #     "HCP-Mn-Ni",
    #     # "HCP-Ni-Fe",
    #     # "HCP-Ni-Mn",
    #     # "HCP-Ni-Ni",
    #     # "SER-Fe",
    #     "SER-Mn",
    #     "SER-Ni",
    # ]
    names = [
        "FCC-Fe-Fe",
    ]
    modeler = Modeler(api_key)

    vasp_args = {
        "handlers": [VaspErrorHandler()],
        "vasp_cmd": "/bin/bash -c 'module load vasp-cpu && srun vasp_std'",
        "vasp_gamma_cmd": "/bin/bash -c 'module load vasp-cpu && srun vasp_gam'",
    }
    incar_settings = {
        "KPAR": 2,
        "NCORE": 1,
        "GGA": "PE",
    }
    potcar_functional = "PBE_64"

    for name in names:
        workdir = Path("data/EndMembers").joinpath(name)
        posdir = Path("data/poscars")
        flowdir = workdir.joinpath("qhaflow")
        json_path = workdir.joinpath(f"{name}-qha.json")
        workdir.mkdir(parents=True, exist_ok=True)

        # # Skip if already done
        # if json_path.exists():
        #     with open(json_path, "r", encoding="utf-8") as jf:
        #         result = json.load(jf)
        #     if result.get("state") == "successful":
        #         log.info(f"System {name} already done")
        #         continue
        
        log.info(f"System {name} started")

        # Get structure
        posdir.mkdir(parents=True, exist_ok=True)
        struct = modeler.get_poscar(name, posdir)
        elems = sorted(str(el) for el in struct.composition.elements)
        log.info(f"System {name} \n sturcture {struct.as_dict()} \n {elems}")

        # Run QHA workflow
        if flowdir.exists():
            shutil.rmtree(flowdir)
        flowdir.mkdir(parents=True, exist_ok=True)

        try:
            worker = QhaWorker(vasp_args, potcar_functional, incar_settings)
            qha_data = worker.run_qha(name, struct, flowdir)
            if not qha_data:
                result = {"name": name, "state": "failed", "struct": struct.as_dict()}
            qha_struct = Structure.from_dict(qha_data.get("structure"))
            qha_elems = sorted(str(el) for el in qha_struct.composition.elements)
            log.info(f"System {name} \n qha structure {qha_struct.as_dict()} \n {qha_elems}")
            if elems != qha_elems:
                result = {
                    "name": name,
                    "state": "failed",
                    "struct": struct.as_dict(),
                    "qha_struct": qha_struct.as_dict(),
                    **qha_data,
                }

            result = {"name": name, "state": "successful", **qha_data}

        except Exception as e:
            log.error(f"System {name} failed: {e}")
            result = {"name": name, "state": "failed"}

        with open(json_path, "w", encoding="utf-8") as jf:
            json.dump(result, jf, indent=2, cls=DateTimeEncoder)

        log.info(f"System {name} done")


if __name__ == "__main__":
    main()

# %%
