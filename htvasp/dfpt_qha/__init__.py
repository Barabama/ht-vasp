"""HT-VASP DFPT-QHA module.

DFPT (``IBRION=8``) phonons, an exact-commensurate-q stability verdict, and a
quasi-harmonic analysis, built on atomate2 makers and jobflow.

Main components:
- DfptQhaMaker: EOS -> per-volume DFPT -> classification -> thermal -> QHA
- DfptVaspMaker / DfptRelaxSetGenerator / DfptInputSetGenerator: VASP inputs
- DfptTaskDoc / DfptForceConstantsProvenance: the collector's output contract
- diagnose_stability / DfptStabilityReport: tier-A checks and tier-B soft modes
- audit_vasp_inputs: what custodian rewrote, original against actual
- analyze_dfpt_qha / dfpt_thermal_properties / volume_free_energy: F(V,T) and G(T)
- force_constants_from_vasprun / force_constants_from_dynmat: VASP's DFPT output

Note: DfptQhaWorker is in htvasp.workflows.dfpt_qha.
"""

from htvasp.dfpt_qha.flow_jobs import (
    classify_dfpt_stability,
    collect_validated_volumes,
    spawn_dfpt_volume_jobs,
)
from htvasp.dfpt_qha.input_audit import audit_vasp_inputs
from htvasp.dfpt_qha.input_set import DfptInputSetGenerator, DfptRelaxSetGenerator
from htvasp.dfpt_qha.jobs import (
    DfptRelaxMaker,
    DfptVaspMaker,
    collect_dfpt_task_doc,
    force_constants_from_dynmat,
    force_constants_from_vasprun,
    parse_dynmat,
    validate_dfpt_task_doc,
)
from htvasp.dfpt_qha.maker import DfptQhaMaker
from htvasp.dfpt_qha.phonons import (
    as_structure,
    assert_same_ordering,
    ordering_fingerprint,
    phonopy_for,
    supercell_structure,
)
from htvasp.dfpt_qha.qha import (
    analyze_dfpt_qha,
    birch_murnaghan_3rd,
    dfpt_thermal_properties,
    volume_free_energy,
)
from htvasp.dfpt_qha.stability import (
    DfptStabilityReport,
    commensurate_qpoints,
    diagnose_stability,
    force_constants_checks,
    force_constants_sha256,
)
from htvasp.dfpt_qha.task_doc import DfptForceConstantsProvenance, DfptTaskDoc

__all__ = [
    "DfptForceConstantsProvenance",
    "DfptInputSetGenerator",
    "DfptQhaMaker",
    "DfptRelaxMaker",
    "DfptRelaxSetGenerator",
    "DfptStabilityReport",
    "DfptTaskDoc",
    "DfptVaspMaker",
    "analyze_dfpt_qha",
    "as_structure",
    "assert_same_ordering",
    "audit_vasp_inputs",
    "birch_murnaghan_3rd",
    "classify_dfpt_stability",
    "collect_dfpt_task_doc",
    "collect_validated_volumes",
    "commensurate_qpoints",
    "dfpt_thermal_properties",
    "diagnose_stability",
    "force_constants_checks",
    "force_constants_from_dynmat",
    "force_constants_from_vasprun",
    "force_constants_sha256",
    "ordering_fingerprint",
    "parse_dynmat",
    "phonopy_for",
    "spawn_dfpt_volume_jobs",
    "supercell_structure",
    "validate_dfpt_task_doc",
    "volume_free_energy",
]
