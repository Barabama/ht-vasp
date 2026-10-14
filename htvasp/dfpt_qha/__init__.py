"""DFPT-QHA endpoint pool: recipes, structures, gates and phonon access.

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
    "read_incar": ("htvasp.dfpt_qha.gate", "read_incar"),
    "read_kpoints": ("htvasp.dfpt_qha.gate", "read_kpoints"),
    "read_poscar": ("htvasp.dfpt_qha.gate", "read_poscar"),
    "read_potcar": ("htvasp.dfpt_qha.gate", "read_potcar"),
    "run_gate": ("htvasp.dfpt_qha.gate", "run_gate"),
    "validate_dir": ("htvasp.dfpt_qha.gate", "validate_dir"),
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
    "DryRunResult",
    "ELEM_MPID_SUPPLEMENT",
    "EndpointRecord",
    "GateResult",
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
    "Smearing",
    "Validation",
    "axis_contains_gamma",
    "build_pool",
    "degeneracy",
    "dryrun_case",
    "eigenvector",
    "enumerate_endpoints",
    "gamma_centred",
    "iter_cases",
    "load_pair",
    "read_incar",
    "read_kpoints",
    "read_poscar",
    "read_potcar",
    "run_gate",
    "validate_dir",
]
