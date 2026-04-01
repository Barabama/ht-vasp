"""
HT-VASP OJ Module

OstravaJ magnetic exchange calculation components.

Main components:
- OJConfig: Configuration (pydantic model)
- OJInputSetGenerator: Input set generator (inherits VaspInputGenerator)
- OJMaker: atomate2-style Maker
- OJResult: Output schema
- Job functions: oj_generate, oj_vasp, oj_solve, oj_workflow
- write_oj_input_set: Input file generation

Note: OJWorker is in htvasp.workflows.oj module
"""

from htvasp.oj.config import OJConfig
from htvasp.oj.input_set import OJInputSetGenerator, write_oj_input_set
from htvasp.oj.jobs import oj_generate, oj_vasp, oj_solve, oj_workflow
from htvasp.oj.maker import OJMaker, OJSimpleMaker, OJRelaxMaker
from htvasp.oj.task_doc import OJResult

__all__ = [
    "OJConfig",
    "OJInputSetGenerator",
    "write_oj_input_set",
    "oj_generate",
    "oj_vasp",
    "oj_solve",
    "oj_workflow",
    "OJMaker",
    "OJSimpleMaker",
    "OJRelaxMaker",
    "OJResult",
]
