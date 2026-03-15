#%%
from atomate2.vasp.jobs.core import RelaxMaker
from atomate2.vasp.sets.core import RelaxSetGenerator
from jobflow.core.store import JobStore
from pymatgen.core import Structure
from maggma.stores import JSONStore
from run_locally_custom import run_locally_custom
import os

vasp_args = {
    "handlers": [],
    "vasp_cmd": "/bin/bash -c 'module load vasp-cpu && mpirun -np 32 vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c 'module load vasp-cpu && mpirun -np 32 vasp_gam'",
}

# %%

store_path = os.path.abspath("test_store.json")
store = JobStore(JSONStore(store_path, read_only=False))
store.connect()

# construct an FCC silicon structure
si_structure = Structure(
    lattice=[[0, 2.73, 2.73], [2.73, 0, 2.73], [2.73, 2.73, 0]],
    species=["Si", "Si"],
    coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
)

# make a relax job to optimise the structure
relax_job = RelaxMaker(
    run_vasp_kwargs=vasp_args,
    input_set_generator=RelaxSetGenerator(
        user_incar_settings={"GGA": "PE"},
    ),
).make(si_structure)

# run the job with custom directory naming
# 使用 job 名作为目录名
run_locally_custom(
    relax_job, 
    store=store, 
    root_dir="test_relax_flow_custom",
    dir_format="{name}",  # 使用 job 名
    dir_prefix="", 
    dir_suffix=""
)

print("\n=== Job completed ===")
print("Check the directory: test_relax_flow_custom/relax/")
