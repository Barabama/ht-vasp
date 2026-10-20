"""Job names shared by the flow builder and the Worker.

The Worker finds the analysis, the stability verdict and a failed validation by
NAME in the job store.  Those names are set in ``flow_jobs`` with f-strings and
queried in ``workflows.dfpt_qha`` with regexes; before this module they were two
independent literals, so renaming a job silently turned ``get_result`` into a
``None``.

They live here, with no imports, so both sides can use them without the Worker
depending on the flow module (or the other way round).
"""

from __future__ import annotations

#: terminal job: EOS, G(T), V(T)
ANALYSIS = "analyze DFPT-QHA"
#: one verdict per validated volume, endpoint gate
STABILITY = "classify DFPT stability"
#: per-volume gate on the collected document (prefix, then "<i>/<n>")
VALIDATION = "validate DFPT volume"
#: per-volume collector (prefix)
COLLECT = "collect DFPT volume"
#: per-volume DFPT VASP job (prefix)
DFPT = "DFPT volume"

__all__ = ["ANALYSIS", "COLLECT", "DFPT", "STABILITY", "VALIDATION"]
