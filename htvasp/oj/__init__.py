"""
HT-VASP OJ Module

OstravaJ magnetic exchange calculation components.

Main components:
- OJConfig: Configuration (pydantic model)
- OJMaker: atomate2-style Maker
- OJResult: Output schema
- Job functions: oj_generate, oj_vasp, oj_solve, oj_workflow
- write_oj_inputs: Input file generation

Note: OJWorker is in htvasp.workflows.oj module
"""

from htvasp.oj.config import OJConfig
from htvasp.oj.input_set import write_oj_inputs
from htvasp.oj.jobs import oj_generate, oj_vasp, oj_solve, oj_workflow
from htvasp.oj.maker import OJMaker, OJSimpleMaker, OJRelaxMaker
from htvasp.oj.task_doc import OJResult

__all__ = [
    "OJConfig",
    "write_oj_inputs",
    "oj_generate",
    "oj_vasp",
    "oj_solve",
    "oj_workflow",
    "OJMaker",
    "OJSimpleMaker",
    "OJRelaxMaker",
    "OJResult",
]
