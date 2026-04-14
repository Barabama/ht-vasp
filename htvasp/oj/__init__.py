"""
HT-VASP OJ Module

OstravaJ magnetic exchange calculation components.

Main components:
- OJInputSetGenerator: Input set generator (inherits VaspInputGenerator)
- OJMaker: atomate2-style Maker with per-flip jobs
- OJResult: Output schema
- Job functions: oj_generate, create_flip_jobs, oj_solve
- write_oj_input_set: Input file generation

Note: OJWorker is in htvasp.workflows.oj module
"""

from htvasp.oj.input_set import OJInputSetGenerator, write_oj_input_set
from htvasp.oj.jobs import oj_generate, oj_solve, create_flip_jobs
from htvasp.oj.maker import OJMaker
from htvasp.oj.task_doc import OJResult

__all__ = [
    "OJInputSetGenerator",
    "write_oj_input_set",
    "oj_generate",
    "oj_solve",
    "create_flip_jobs",
    "OJMaker",
    "OJResult",
]
