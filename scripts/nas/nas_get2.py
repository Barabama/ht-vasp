#!/usr/bin/env python
"""Read-only fetch of NAS files over SFTP using .env credentials.

Usage:
  nas_get2.py cat  <win-path> [<win-path> ...]
  nas_get2.py pull <win-file> <local-file>
Never writes to the NAS.
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


def to_sftp_path(win: str) -> str:
    p = win.replace("\\", "/")
    if len(p) > 1 and p[1] == ":":
        p = "/" + p
    return p


def connect():
    e = load_env()
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(e["NAS_HOST"], port=int(e.get("NAS_SSH_PORT", 22)),
              username=e["NAS_USER"], password=e["NAS_PASS"],
              timeout=25, allow_agent=False, look_for_keys=False)
    return c, c.open_sftp()


def main():
    mode = sys.argv[1]
    c, sftp = connect()
    try:
        if mode == "cat":
            for win in sys.argv[2:]:
                print(f"########## {win}")
                try:
                    with sftp.open(to_sftp_path(win), "rb") as f:
                        data = f.read()
                    try:
                        print(data.decode("utf-8").rstrip())
                    except UnicodeDecodeError:
                        print(f"[binary, {len(data)} bytes]")
                except IOError as ex:
                    print(f"<ERROR: {ex}>")
                print()
        elif mode == "pull":
            src, dst = to_sftp_path(sys.argv[2]), sys.argv[3]
            os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)
            sftp.get(src, dst)
            print(f"{os.path.getsize(dst):>12,}  {dst}")
        else:
            raise SystemExit("mode must be cat|pull")
    finally:
        sftp.close()
        c.close()


if __name__ == "__main__":
    main()
