#!/usr/bin/env python
"""Read-only probe of the MCMF-NAS over SSH, with credentials from the repo .env.

Usage: nas_probe2.py "<cmd>" ["<cmd>" ...]
Never writes, never deletes.
"""
import os
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


def connect():
    e = load_env()
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(e["NAS_HOST"], port=int(e.get("NAS_SSH_PORT", 22)),
              username=e["NAS_USER"], password=e["NAS_PASS"],
              timeout=25, allow_agent=False, look_for_keys=False)
    return c, e


def main():
    c, e = connect()
    cmds = sys.argv[1:] or ["ver"]
    for cmd in cmds:
        _, out, err = c.exec_command(cmd, timeout=180)
        o = out.read().decode(errors="replace")
        er = err.read().decode(errors="replace")
        print(f"$ {cmd}")
        if o:
            print(o.rstrip())
        if er:
            print("[stderr]", er.rstrip())
        print()
    c.close()


if __name__ == "__main__":
    main()
