"""
HT-VASP TB2J Module

TB2J magnetic exchange calculation components (collinear Wannier90-based).

Main components:
- Tb2jInputSetGenerator: Input set generator (inherits VaspInputGenerator),
  writes wannier90.win + Wannier90-interface INCAR keys
- Tb2jMaker: atomate2-style Maker composing W90 SCF -> TB2J solve
- Tb2jResult: Output schema (parsed J pairs, tensors, summary)
- tb2j_solve: Job function running ``python -m TB2J.scripts.wann2J``

Note: the full end-to-end Tb2jWorker (R3 relax -> W90 SCF -> TB2J solve) lives
in ``htvasp.workflows.tb2j``; this subpackage provides the reusable pieces.
"""

from htvasp.tb2j.input_set import Tb2jInputSetGenerator, ORBITALS_PER_ATOM
from htvasp.tb2j.jobs import tb2j_solve
from htvasp.tb2j.maker import Tb2jMaker
from htvasp.tb2j.task_doc import Tb2jResult

__all__ = [
    "Tb2jInputSetGenerator",
    "ORBITALS_PER_ATOM",
    "tb2j_solve",
    "Tb2jMaker",
    "Tb2jResult",
]
