from htvasp.model.sqs import (
    DEFAULT_SUPERCELLS,
    SqsCalculationSpec,
    SqsCompositionPoint,
    SqsPlan,
    SqsStructureArtifact,
    SqsSystem,
    default_supercell,
    payload_checksum,
)
from htvasp.model.templates import (
    TemplateError,
    TemplateRecord,
    TemplateRegistry,
    TemplateSnapshot,
    file_checksum,
    structure_checksum,
)
from htvasp.model.dfpt_qha_recipe import Recipe, Validation, MissingRecipeField, RecipeInconsistent, load_pair
from htvasp.model.dfpt_qha_pool import PoolEndpoint, PoolManifest, POOL_SCHEMA
from htvasp.model.sqs_schema import (
    CalculationProvenance,
    CalphadMapping,
    DatasetValidation,
    EV_PER_ATOM_TO_J_PER_MOL,
    PureReference,
    SCHEMA_VERSION,
    SqsMixingData,
    SqsMixingError,
    SqsMixingPoint,
    SqsProvenance,
    StaticObservation,
    build_bcc_a2_mapping,
    build_sqs_mixing_data,
)


def __getattr__(name: str):
    """Keep the historical Endmember export lazy and network-free at import."""

    if name == "Endmember":
        from htvasp.model.endmember import Endmember

        return Endmember
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "Endmember",
    "DEFAULT_SUPERCELLS",
    "SqsCalculationSpec",
    "SqsCompositionPoint",
    "SqsPlan",
    "SqsStructureArtifact",
    "SqsSystem",
    "default_supercell",
    "payload_checksum",
    "TemplateError",
    "TemplateRecord",
    "TemplateRegistry",
    "TemplateSnapshot",
    "file_checksum",
    "structure_checksum",
    "CalculationProvenance",
    "CalphadMapping",
    "DatasetValidation",
    "EV_PER_ATOM_TO_J_PER_MOL",
    "PureReference",
    "SCHEMA_VERSION",
    "SqsMixingData",
    "SqsMixingError",
    "SqsMixingPoint",
    "SqsProvenance",
    "StaticObservation",
    "build_bcc_a2_mapping",
    "build_sqs_mixing_data",
    "Recipe",
    "Validation",
    "MissingRecipeField",
    "RecipeInconsistent",
    "load_pair",
    "PoolEndpoint",
    "PoolManifest",
    "POOL_SCHEMA",
]
