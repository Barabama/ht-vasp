"""Write one line per completed endpoint: name|cluster|status|applicability."""
import json, glob, sys
from pathlib import Path

cluster = sys.argv[1]
out = sys.argv[2]
rows = set()
for f in glob.glob("data/endpoints/*/dfpt-qhaflow/runs/*/*-dfpt-qha.json"):
    d = json.load(open(f))
    if d.get("ok") is not True:
        continue
    name = d.get("endpoint") or Path(f).name.split("-dfpt-qha")[0]
    rows.add("%s|%s|done %s|%s" % (name, cluster, cluster,
                                   d.get("qha_applicability") or "unknown"))
Path(out).write_text("\n".join(sorted(rows)) + "\n")
print(len(rows))
