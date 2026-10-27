#!/usr/bin/env python
"""Fast structural census of the 16symbols tree: what does each endpoint carry?

Only lists directories, never walks into per-volume trees, so it is quick.
Writes nothing to the NAS.
"""
import json
import os
import posixpath
import stat
import sys

import paramiko

ENV = "/nfs_hdd/2025/gaominliang/.env"


def load_env(path=ENV):
    d = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip().strip('"').strip("'")
    return d


def to_sftp_path(win):
    p = win.replace("\\", "/")
    if len(p) > 1 and p[1] == ":":
        p = "/" + p
    return p


def main():
    e = load_env()
    root_win = sys.argv[1]
    out = sys.argv[2]
    root = to_sftp_path(root_win)

    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(e["NAS_HOST"], port=int(e.get("NAS_SSH_PORT", 22)),
              username=e["NAS_USER"], password=e["NAS_PASS"],
              timeout=25, allow_agent=False, look_for_keys=False)
    sftp = c.open_sftp()

    eps = sorted(a.filename for a in sftp.listdir_attr(root)
                 if stat.S_ISDIR(a.st_mode))
    print(f"{len(eps)} endpoint directories")
    rows = []
    for i, ep in enumerate(eps, 1):
        p = posixpath.join(root, ep)
        try:
            entries = sftp.listdir_attr(p)
        except IOError as ex:
            rows.append({"endpoint": ep, "error": str(ex)})
            continue
        names = {a.filename for a in entries}
        dirs = {a.filename for a in entries if stat.S_ISDIR(a.st_mode)}
        files = names - dirs
        qha_dirs = sorted(d for d in dirs if d.startswith("QHA-"))
        dfpt_vols = []
        dfpt_has = "dfpt" in dirs
        if dfpt_has:
            try:
                sub = sftp.listdir_attr(posixpath.join(p, "dfpt"))
                dfpt_vols = sorted(a.filename for a in sub
                                   if stat.S_ISDIR(a.st_mode) and a.filename.startswith("QHA-"))
            except IOError:
                pass
        r3_has = "R3" in dirs
        r3_files = []
        if r3_has:
            try:
                r3_files = sorted(a.filename for a in sftp.listdir_attr(posixpath.join(p, "R3")))
            except IOError:
                pass
        tp_top = sorted(f for f in files if f.startswith("thermal_properties-"))
        rows.append({
            "endpoint": ep,
            "has_dfpt": dfpt_has,
            "dfpt_volumes": dfpt_vols,
            "has_R3": r3_has,
            "R3_files": r3_files,
            "curve_dirs": sorted(d for d in dirs if not d.startswith("QHA-")),
            "top_files": sorted(files),
            "top_thermal_properties": tp_top,
            "qha_dirs": qha_dirs,
        })
        if i % 40 == 0:
            print(f"  {i}/{len(eps)}", flush=True)

    with open(out, "w", encoding="utf-8") as f:
        json.dump({"schema": "nas16-structure-census/v1", "root": root_win,
                   "n_endpoints": len(rows), "endpoints": rows},
                  f, indent=1)
    print("->", out)
    sftp.close()
    c.close()


if __name__ == "__main__":
    main()
