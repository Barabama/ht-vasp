#!/usr/bin/env python
"""Small read/write SSH+SFTP helper for the MCMF-NAS, credentials from .env.

Writes are refused unless the path is inside NAS_SHARED_ENDMEMBERS_PATH
(the one tree the user authorised). Everything else stays read-only.

Usage:
  nas.py run  "<cmd>" ["<cmd>" ...]
  nas.py put  <local> <win-remote>
  nas.py get  <win-remote> <local>
  nas.py ls   <win-remote>
  nas.py mkdir <win-remote>
"""
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


E = load_env()
AUTH_ROOT = E["NAS_SHARED_ENDMEMBERS_PATH"]          # \e\nas-shared\Endmembers\16symbols-cxy-hmd


def to_sftp_path(win: str) -> str:
    """Windows path -> SFTP path. Accepts both `E:\\a\\b` and `\\e\\a\\b`."""
    p = win.replace("\\", "/")
    if len(p) > 1 and p[1] == ":":
        p = "/" + p
    return p


def _canon(win: str) -> str:
    """Drive letter with or without its colon, lowercased, for comparison."""
    p = posixpath.normpath(to_sftp_path(win))
    if len(p) > 2 and p[0] == "/" and p[2] == ":":
        p = "/" + p[1] + p[3:]
    return p.lower()


def sftp_norm(win: str) -> str:
    """Normalised sftp path for the authorised tree."""
    return _canon(AUTH_ROOT)


def assert_writable(win: str) -> None:
    p = _canon(win)
    root = sftp_norm(AUTH_ROOT)
    if not (p == root or p.startswith(root + "/")):
        raise SystemExit(
            f"REFUSED: {win}\n  is outside NAS_SHARED_ENDMEMBERS_PATH ({AUTH_ROOT}); "
            f"this session may only write inside that tree")


def connect():
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(E["NAS_HOST"], port=int(E.get("NAS_SSH_PORT", 22)),
              username=E["NAS_USER"], password=E["NAS_PASS"],
              timeout=25, allow_agent=False, look_for_keys=False)
    return c, c.open_sftp()


def run(c, cmd, timeout=1800):
    _, out, err = c.exec_command(cmd, timeout=timeout)
    o = out.read().decode(errors="replace")
    er = err.read().decode(errors="replace")
    return o, er


def main():
    mode = sys.argv[1]
    c, sftp = connect()
    try:
        if mode == "run":
            for cmd in sys.argv[2:]:
                o, er = run(c, cmd)
                print(f"$ {cmd}")
                if o:
                    print(o.rstrip())
                if er:
                    print("[stderr]", er.rstrip())
                print()
        elif mode == "put":
            local, remote = sys.argv[2], sys.argv[3]
            assert_writable(remote)
            rp = to_sftp_path(remote)
            d = posixpath.dirname(rp)
            try:
                sftp.stat(d)
            except IOError:
                sftp.mkdir(d)
            sftp.put(local, rp)
            print(f"{os.path.getsize(local):>12,}  -> {remote}")
        elif mode == "get":
            remote, local = sys.argv[2], sys.argv[3]
            os.makedirs(os.path.dirname(os.path.abspath(local)), exist_ok=True)
            sftp.get(to_sftp_path(remote), local)
            print(f"{os.path.getsize(local):>12,}  <- {remote}")
        elif mode == "ls":
            for e in sorted(sftp.listdir_attr(to_sftp_path(sys.argv[2])),
                            key=lambda a: a.filename):
                kind = "DIR " if stat.S_ISDIR(e.st_mode) else "    "
                print(f"{kind}{e.st_size:>13,}  {e.filename}")
        elif mode == "mkdir":
            target = sys.argv[2]
            assert_writable(target)
            rp = to_sftp_path(target)
            parts, cur = rp.strip("/").split("/"), ""
            for p in parts:
                cur += "/" + p
                try:
                    sftp.stat(cur)
                except IOError:
                    sftp.mkdir(cur)
                    print("created", cur)
        else:
            raise SystemExit("mode must be run|put|get|ls|mkdir")
    finally:
        sftp.close()
        c.close()


if __name__ == "__main__":
    main()
