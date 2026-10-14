"""DFPT-QHA endpoint pool: recipes, structures, gates, inputs and phonons.

The public surface is re-exported lazily so that importing this package does
not pull pymatgen or phonopy into a process that only wants the recipe types.
"""

from importlib import import_module

_EXPORTS = {
    # recipe contract
    "EndpointRecord": ("htvasp.dfpt_qha.recipe", "EndpointRecord"),
    "Kpoints": ("htvasp.dfpt_qha.recipe", "Kpoints"),
    "MissingRecipeField": ("htvasp.dfpt_qha.recipe", "MissingRecipeField"),
    "Potcar": ("htvasp.dfpt_qha.recipe", "Potcar"),
    "Recipe": ("htvasp.dfpt_qha.recipe", "Recipe"),
    "RecipeError": ("htvasp.dfpt_qha.recipe", "RecipeError"),
    "RecipeInconsistent": ("htvasp.dfpt_qha.recipe", "RecipeInconsistent"),
    "SCHEMA_VERSION": ("htvasp.dfpt_qha.recipe", "SCHEMA_VERSION"),
    "Smearing": ("htvasp.dfpt_qha.recipe", "Smearing"),
    "Validation": ("htvasp.dfpt_qha.recipe", "Validation"),
    "axis_contains_gamma": ("htvasp.dfpt_qha.recipe", "axis_contains_gamma"),
    "gamma_centred": ("htvasp.dfpt_qha.recipe", "gamma_centred"),
    "load_pair": ("htvasp.dfpt_qha.recipe", "load_pair"),
    # endpoint pool
    "ELEM_MPID_SUPPLEMENT": ("htvasp.dfpt_qha.pool", "ELEM_MPID_SUPPLEMENT"),
    "POOL_ELEMENTS": ("htvasp.dfpt_qha.pool", "POOL_ELEMENTS"),
    "PoolEndpoint": ("htvasp.dfpt_qha.pool", "PoolEndpoint"),
    "PoolManifest": ("htvasp.dfpt_qha.pool", "PoolManifest"),
    "build_pool": ("htvasp.dfpt_qha.pool", "build_pool"),
    "enumerate_endpoints": ("htvasp.dfpt_qha.pool", "enumerate_endpoints"),
    # VASP input gate
    "DryRunResult": ("htvasp.dfpt_qha.gate", "DryRunResult"),
    "GateResult": ("htvasp.dfpt_qha.gate", "GateResult"),
    "dryrun_case": ("htvasp.dfpt_qha.gate", "dryrun_case"),
    "iter_cases": ("htvasp.dfpt_qha.gate", "iter_cases"),
    "parse_repeated_floats": ("htvasp.dfpt_qha.gate", "parse_repeated_floats"),
    "read_incar": ("htvasp.dfpt_qha.gate", "read_incar"),
    "read_kpoints": ("htvasp.dfpt_qha.gate", "read_kpoints"),
    "read_poscar": ("htvasp.dfpt_qha.gate", "read_poscar"),
    "read_potcar": ("htvasp.dfpt_qha.gate", "read_potcar"),
    "run_gate": ("htvasp.dfpt_qha.gate", "run_gate"),
    "validate_dir": ("htvasp.dfpt_qha.gate", "validate_dir"),
    # input generation
    "DFPT_ELECTRONIC": ("htvasp.dfpt_qha.input_set", "DFPT_ELECTRONIC"),
    "InputSetError": ("htvasp.dfpt_qha.input_set", "InputSetError"),
    "SHARED_INCAR": ("htvasp.dfpt_qha.input_set", "SHARED_INCAR"),
    "STAGES": ("htvasp.dfpt_qha.input_set", "STAGES"),
    "STATIC_ELECTRONIC": ("htvasp.dfpt_qha.input_set", "STATIC_ELECTRONIC"),
    "incar_for": ("htvasp.dfpt_qha.input_set", "incar_for"),
    "kpoints_for": ("htvasp.dfpt_qha.input_set", "kpoints_for"),
    "magmom_for": ("htvasp.dfpt_qha.input_set", "magmom_for"),
    "mesh_for": ("htvasp.dfpt_qha.input_set", "mesh_for"),
    "potcar_symbols": ("htvasp.dfpt_qha.input_set", "potcar_symbols"),
    "round_even": ("htvasp.dfpt_qha.input_set", "round_even"),
    "supercell_kpoints_for": ("htvasp.dfpt_qha.input_set", "supercell_kpoints_for"),
    "write_stage_inputs": ("htvasp.dfpt_qha.input_set", "write_stage_inputs"),
    # phonons
    "BandEigenvector": ("htvasp.dfpt_qha.phonons", "BandEigenvector"),
    "degeneracy": ("htvasp.dfpt_qha.phonons", "degeneracy"),
    "eigenvector": ("htvasp.dfpt_qha.phonons", "eigenvector"),
}


def __getattr__(name: str):
    try:
        module_name, attribute = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc
    value = getattr(import_module(module_name), attribute)
    globals()[name] = value
    return value


__all__ = [
    "BandEigenvector",
    "DFPT_ELECTRONIC",
    "DryRunResult",
    "ELEM_MPID_SUPPLEMENT",
    "EndpointRecord",
    "GateResult",
    "InputSetError",
    "Kpoints",
    "MissingRecipeField",
    "POOL_ELEMENTS",
    "PoolEndpoint",
    "PoolManifest",
    "Potcar",
    "Recipe",
    "RecipeError",
    "RecipeInconsistent",
    "SCHEMA_VERSION",
    "SHARED_INCAR",
    "STAGES",
    "STATIC_ELECTRONIC",
    "Smearing",
    "Validation",
    "axis_contains_gamma",
    "build_pool",
    "degeneracy",
    "dryrun_case",
    "eigenvector",
    "enumerate_endpoints",
    "gamma_centred",
    "incar_for",
    "iter_cases",
    "kpoints_for",
    "load_pair",
    "magmom_for",
    "mesh_for",
    "parse_repeated_floats",
    "potcar_symbols",
    "read_incar",
    "read_kpoints",
    "read_poscar",
    "read_potcar",
    "round_even",
    "run_gate",
    "supercell_kpoints_for",
    "validate_dir",
    "write_stage_inputs",
]
