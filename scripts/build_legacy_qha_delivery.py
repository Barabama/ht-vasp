#!/usr/bin/env python
"""Materialise an A-line -> B-line handoff tree for an ingested legacy QHA endpoint.

READ-ONLY with respect to the source. Nothing under ``legacy-as-received/`` is
written, moved or deleted: every file in the delivery is a byte-exact copy, and
each copy's sha256 is checked against the mirror's own manifest before the
delivery is declared complete.

What the delivery is
--------------------
``em-tdb-v2``'s plan (§B.5.7) explicitly refuses to freeze a traditional-DFPT
directory layout until the user has the real data in hand, and (§B.5.8) refuses
to let a directory name imply physical semantics. So this tool does **not**
invent a canonical layout and call it a contract. It produces:

* an A-line *proposal* for the layout, with every choice tied to the evidence
  that supports it and to the ``§B.5.7`` checklist item it answers;
* the observables themselves, copied verbatim, so the B-line has real bytes;
* the fail-closed audit report, so the B-line can compute an applicability
  status without re-deriving anything;
* pointers (path + sha256) into the verbatim mirror for the raw VASP layer,
  which is deliberately not duplicated.

Usage:
    build_legacy_qha_delivery.py --root data/legacy-dfpt-qha-pilot/BCC-Co-Ti \\
        --outdir data/legacy-dfpt-qha-delivery [--label BCC-Co-Ti-2]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
from legacy_qha_audit import read_v_e, sha256_file  # noqa: E402

# The raw VASP layer: referenced by hash, never duplicated into the delivery.
RAW_GLOBS = ("FORCE_CONSTANTS", "OUTCAR", "vasprun.xml", "log-DFPT", "INCAR",
             "INCAR-R3f", "OSZICAR", "DYNMAT", "disp.yaml", "phonopy.yaml",
             "phonopy_disp.yaml", "mesh.yaml", "POSCAR", "POSCAR-unitcell",
             "CONTCAR", "REPORT")


def copy_verified(src: Path, dst: Path, manifest_hashes: dict[str, str],
                  rel_src: str) -> dict[str, Any]:
    """Copy one file byte-exactly and prove it against the mirror manifest."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        raise SystemExit(
            f"refusing to overwrite an existing delivery file: {dst}\n"
            f"(the delivery tree is meant to be created once; move it aside "
            f"if you really want to regenerate it)")
    shutil.copy2(src, dst)
    got = sha256_file(dst)
    if src.stat().st_size != dst.stat().st_size:
        raise SystemExit(f"size mismatch copying {src} -> {dst}")
    expect = manifest_hashes.get(rel_src)
    return {
        "delivery_path": str(dst),
        "mirror_path": rel_src,
        "size_bytes": dst.stat().st_size,
        "sha256": got,
        "matches_mirror_manifest": (None if expect is None else got == expect),
    }


def load_manifest_hashes(root: Path) -> dict[str, str]:
    mp = root / "source-manifest.json"
    if not mp.is_file():
        return {}
    d = json.loads(mp.read_text())
    return {f["relative_path"]: f["sha256"] for f in d.get("files", [])}


def raw_index(rec: Path, subdirs: list[Path]) -> list[dict[str, Any]]:
    """Path+hash pointers into the mirror's VASP layer. No bytes are copied."""
    out = []
    for d in subdirs:
        if not d.is_dir():
            continue
        for name in RAW_GLOBS:
            p = d / name
            if p.is_file():
                out.append({
                    "mirror_path": str(p.relative_to(rec)),
                    "size_bytes": p.stat().st_size,
                    "sha256": sha256_file(p),
                })
    return sorted(out, key=lambda r: r["mirror_path"])


def temperature_grid(path: Path) -> list[float] | None:
    """[T_first, T_last, T_step, n_points] for a two-column numeric curve."""
    rows = []
    for line in path.read_text().splitlines():
        s = line.split()
        if not s or s[0].startswith("#"):
            continue
        try:
            rows.append((float(s[0]), float(s[1])))
        except (ValueError, IndexError):
            continue
    if len(rows) < 2:
        return None
    return [rows[0][0], rows[-1][0], rows[1][0] - rows[0][0], len(rows)]


def build(root: Path, outdir: Path, label: str,
          external_structure_json: Path | None = None,
          external_structure_template: Path | None = None) -> dict[str, Any]:
    rec = root / "legacy-as-received"
    if not rec.is_dir():
        raise SystemExit(f"missing {rec}")
    audit_path = root / "audit" / "audit.json"
    if not audit_path.is_file():
        raise SystemExit(
            f"missing {audit_path}; run scripts/legacy_qha_audit.py first")
    audit = json.loads(audit_path.read_text())
    mh = load_manifest_hashes(root)

    ep = outdir / label
    copied: list[dict[str, Any]] = []

    # ---- curves -------------------------------------------------------
    # Some mirrors keep the curves at the root, some inside QHA-<name>/, some
    # in both. Flatten to curves/<basename>, keep the root copy when both
    # exist, and record the alternate mirror path so nothing is silently lost.
    curve_dir = ep / "curves"
    curve_srcs: dict[str, list[Path]] = {}
    for src in sorted(rec.glob("*.dat")) + sorted(rec.glob("QHA-*/*.dat")):
        curve_srcs.setdefault(src.name, []).append(src)
    curve_dupes: dict[str, list[dict[str, Any]]] = {}
    for name, srcs in sorted(curve_srcs.items()):
        root_srcs = [s for s in srcs if s.parent == rec]
        primary = root_srcs[0] if root_srcs else srcs[0]
        copied.append(copy_verified(primary, curve_dir / name, mh,
                                    str(primary.relative_to(rec))))
        others = [s for s in srcs if s != primary]
        if others:
            curve_dupes[name] = [{
                "mirror_path": str(s.relative_to(rec)),
                "sha256": sha256_file(s),
                "identical_to_copied": sha256_file(s) == sha256_file(primary),
            } for s in others]

    # ---- phonon layer -------------------------------------------------
    phon_dir = ep / "phonons"
    phonon_yaml: list[dict[str, Any]] = []
    for src in sorted(rec.glob("thermal_properties-*.yaml")):
        rel = str(src.relative_to(rec))
        copied.append(copy_verified(src, phon_dir / src.name, mh, rel))
        phonon_yaml.append({"file": src.name, "mirror_path": rel})
    for src in sorted(rec.glob("dfpt/QHA-*/thermal_properties-*.yaml")):
        rel = str(src.relative_to(rec))
        tag = src.parent.name.replace("QHA-", "")
        name = f"thermal_properties-from-dfpt-{tag}.yaml"
        copied.append(copy_verified(src, phon_dir / name, mh, rel))
        phonon_yaml.append({"file": name, "mirror_path": rel})

    # ---- structure ----------------------------------------------------
    struct_dir = ep / "structure"
    struct_notes: list[str] = []
    anchor = rec / "POSCAR-unitcell"
    if anchor.is_file():
        copied.append(copy_verified(anchor, struct_dir / "POSCAR-unitcell",
                                    mh, "POSCAR-unitcell"))
        struct_notes.append(
            "POSCAR-unitcell copied from the legacy mirror (shipped at its root)")
    # The relaxed anchor cell, when the mirror kept the R3 relaxation.
    contcar = rec / "R3" / "CONTCAR"
    if contcar.is_file():
        copied.append(copy_verified(contcar, struct_dir / "relaxed-CONTCAR",
                                    mh, "R3/CONTCAR"))
        struct_notes.append(
            "R3/CONTCAR copied: this is the ISIF=3-relaxed cell the volume grid "
            "is built around")
    # The mirror's own root-level cell files, if it has any.
    for name in ("POSCAR", "POSCAR-org", "POSCAR-unitcell"):
        src = rec / name
        if src.is_file():
            copied.append(copy_verified(src, struct_dir / name, mh, name))
            struct_notes.append(f"{name} copied from the legacy mirror root")
    per_vol = sorted(rec.glob("dfpt/QHA-*/POSCAR-unitcell"))
    for src in per_vol:
        rel = str(src.relative_to(rec))
        name = f"volumes/POSCAR-unitcell-QHA-{src.parent.name.replace('QHA-', '')}"
        copied.append(copy_verified(src, struct_dir / name, mh, rel))
    if per_vol:
        struct_notes.append(
            f"{len(per_vol)} per-volume POSCAR-unitcell files copied; these are "
            f"the exact cells the DFPT force constants belong to")
        if not anchor.is_file():
            struct_notes.append(
                "the mirror has NO POSCAR-unitcell at its root; the anchor cell "
                "is therefore only available as the per-volume files above plus "
                "R3/CONTCAR")
    if not anchor.is_file() and not per_vol:
        struct_notes.append(
            "NO structure file anywhere in the legacy mirror: this endpoint "
            "carries the curves and the thermal properties only, so any cell "
            "used to read them has to be supplied from outside")

    # ---- structure supplied from outside the legacy mirror -------------
    if external_structure_json is not None:
        from pymatgen.core import Structure
        d = json.loads(Path(external_structure_json).read_text())
        st = Structure.from_dict(d["structure"])
        ext = struct_dir / "from-ht-vasp"
        ext.mkdir(parents=True, exist_ok=True)
        poscar = ext / "POSCAR-unitcell"
        if poscar.exists():
            raise SystemExit(f"refusing to overwrite {poscar}")
        st.to(filename=str(poscar), fmt="poscar")
        struct_notes.append(
            f"structure/from-ht-vasp/POSCAR-unitcell was generated from "
            f"{external_structure_json} (ht-vasp relaxed cell, "
            f"a={st.lattice.a:.6f} A, V={st.volume:.6f} A^3); it is NOT from the "
            f"legacy mirror and carries no legacy provenance")
        if external_structure_template is not None:
            tsrc = Path(external_structure_template)
            copied.append(copy_verified(tsrc, ext / tsrc.name, {},
                                        str(tsrc)))
            struct_notes.append(
                f"structure/from-ht-vasp/{tsrc.name} is the ht-vasp template "
                f"POSCAR; its lattice constant is a shared placeholder used by "
                f"every BCC-* template in that directory and must not be read "
                f"as this endpoint's cell")

    # ---- recipe layer (small text inputs) ------------------------------
    # The DFPT recipe is what §B.5.7 item 9 asks about ("how is this route's
    # producer / workflow identity proven"), so it is copied rather than
    # merely pointed at. POTCAR is not copied: it is large and licensed, and
    # the B-line needs its identity, not its bytes.
    recipe_dir = ep / "recipe"
    recipe_copied: list[str] = []
    for pat in ("INCAR*", "KPOINTS*", "mesh.conf"):
        for src in sorted(rec.glob(pat)):
            if not src.is_file():
                continue
            copied.append(copy_verified(src, recipe_dir / src.name, mh,
                                        str(src.relative_to(rec))))
            recipe_copied.append(src.name)
    potcar = rec / "POTCAR"
    if potcar.is_file():
        lines = [l.rstrip() for l in potcar.read_text(errors="replace").splitlines()
                 if re.match(r"\s*(TITEL|POMASS|ZVAL|VRHFIN)", l)]
        (recipe_dir / "POTCAR_IDENTITY.txt").write_text(
            "identity of the POTCAR blocks, in file order.\n"
            "The full POTCAR is not copied; it stays in the mirror at "
            "legacy-as-received/POTCAR.\n"
            "Blocks beyond the number of species POSCAR declares are never "
            "read by VASP and are template residue.\n\n"
            + "\n".join(lines) + "\n", encoding="utf-8")
        recipe_copied.append("POTCAR_IDENTITY.txt")

    # ---- audit + raw index --------------------------------------------
    copied.append(copy_verified(audit_path, ep / "audit" / "audit.json", {}, "audit/audit.json"))

    subdirs = [rec / "R3"] + sorted(rec.glob("dfpt/QHA-*")) + \
              sorted(rec.glob("dfpt/QHA-*/R3f"))
    raw = raw_index(rec, subdirs)
    for src in sorted(rec.glob("POTCAR*")):
        if src.is_file():
            raw.append({"mirror_path": src.name, "size_bytes": src.stat().st_size,
                        "sha256": sha256_file(src)})
    raw.sort(key=lambda r: r["mirror_path"])
    raw_dir = ep / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    (raw_dir / "RAW_INDEX.json").write_text(json.dumps({
        "schema": "legacy-qha-raw-index/v1",
        "note": ("pointers into legacy-as-received/, by path and sha256. The raw "
                 "VASP layer is deliberately not duplicated into the delivery; "
                 "these are the files that would have to be re-read to "
                 "re-derive the force constants."),
        "mirror_root": str(rec.resolve()),
        "n_files": len(raw),
        "total_bytes": sum(r["size_bytes"] for r in raw),
        "files": raw,
    }, indent=2) + "\n", encoding="utf-8")

    # ---- curve metadata ----------------------------------------------
    v_e = rec / "v-e.dat"
    V, E = read_v_e(v_e)
    curves = {
        "v_e": {
            "file": "curves/v-e.dat",
            "sha256": sha256_file(v_e),
            "columns": ["volume_angstrom3", "energy_eV"],
            "n_points": int(len(V)),
            "volume_range_angstrom3": [float(V.min()), float(V.max())],
            "volume_step_angstrom3": float(V[1] - V[0]),
            "energy_range_eV": [float(E.min()), float(E.max())],
            "normalization": audit["eos_input"]["normalization"],
            "provenance": audit["eos_input"]["provenance"],
        },
    }
    for name in ("gibbs-temperature.dat", "helmholtz-volume.dat", "Cv-volume.dat",
                 "entropy-volume.dat", "Cp-temperature.dat",
                 "Cp-temperature_polyfit.dat", "bulk_modulus-temperature.dat",
                 "dsdv-temperature.dat", "gruneisen-temperature.dat",
                 "thermal_expansion.dat", "volume-temperature.dat",
                 "volume_expansion.dat"):
        for cand in (rec / name, rec / "QHA-TiCo" / name):
            if cand.is_file():
                g = temperature_grid(cand)
                if g:
                    curves[name] = {
                        "file": "curves/" + name,
                        "sha256": sha256_file(cand),
                        "n_points": g[3],
                        "first_temperature_K": g[0],
                        "last_temperature_K": g[1],
                        "temperature_step_K": g[2],
                    }
                break

    summary = {
        "schema": "legacy-qha-delivery-inventory/v1",
        "generated_by": "scripts/build_legacy_qha_delivery.py",
        "label": label,
        "endpoint_root": str(root.resolve()),
        "mirror_root": str(rec.resolve()),
        "audit": {
            "status": audit["status"],
            "vasp_layer_present": audit.get("vasp_layer_present"),
            "n_failures": len(audit["failures"]),
            "n_warnings": len(audit.get("warnings", [])),
            "failures": audit["failures"],
            "warnings": audit.get("warnings", []),
            "notes": audit.get("notes", []),
            "report": "audit/audit.json",
        },
        "curves": curves,
        "curve_alternate_copies": curve_dupes,
        "phonon_yaml_copied": phonon_yaml,
        "recipe_files_copied": recipe_copied,
        "structure_notes": struct_notes,
        "raw_index": {"file": "raw/RAW_INDEX.json", "n_files": len(raw),
                      "total_bytes": sum(r["size_bytes"] for r in raw)},
        "files": copied,
    }
    # Record every copied file's location relative to the delivery root, so the
    # inventory stays valid wherever the tree is moved to.
    for f in summary["files"]:
        f["delivery_path"] = str(
            Path(f["delivery_path"]).resolve().relative_to(outdir.resolve()))
    return summary

def overlay_content(content_dir: Path, ep: Path, label: str,
                    generated: set[Path]) -> list[str]:
    """Copy the hand-authored documents in, so a rebuild regenerates everything.

    Shared files live at the content root (README.md, LAYOUT_PROPOSAL.md); the
    rest live under ``<content_dir>/<label>/``. Nothing here is generated, so
    the overlay may refresh its own files but must never land on a path the
    mechanical part of this tool produced.
    """
    written: list[str] = []
    if not content_dir.is_dir():
        return written
    pairs: list[tuple[Path, Path]] = []
    for src in sorted(content_dir.iterdir()):
        if src.is_file():
            pairs.append((src, ep.parent / src.name))
    per = content_dir / label
    if per.is_dir():
        for src in sorted(per.rglob("*")):
            if src.is_file():
                pairs.append((src, ep / src.relative_to(per)))
    for src, dst in pairs:
        if dst in generated:
            raise SystemExit(
                f"authored content {src} would overwrite a generated file {dst}")
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        written.append(str(dst.relative_to(ep.parent)))
    return written


def build_index(outdir: Path) -> dict[str, Any]:
    """Root-level DELIVERY.json, assembled from what is actually on disk."""
    endpoints = {}
    for ep in sorted(p for p in outdir.iterdir() if p.is_dir() and p.name[0] != "_"):
        dj = ep / "DATASET.json"
        iv = ep / "INVENTORY.json"
        if not (dj.is_file() and iv.is_file()):
            continue
        d = json.loads(dj.read_text())
        inv = json.loads(iv.read_text())
        endpoints[ep.name] = {
            "nas_path": d["label_provenance"]["nas_path"],
            "endpoint_content_digest_sha256": d["source"]["endpoint_content_digest_sha256"],
            "audit_status": d["audit"]["status"],
            "audit_n_failures": d["audit"]["n_failures"],
            "audit_n_warnings": d["audit"]["n_warnings"],
            "applicability": d["claims"]["applicability"]["value"],
            "vasp_layer_present": inv["audit"]["vasp_layer_present"],
            "mirror_root": d["source"]["mirror_root"],
            "copied_files": len(inv["files"]),
            "raw_vasp_files_referenced": inv["raw_index"]["n_files"],
            "claim_status": {k: v.get("status") for k, v in d["claims"].items()},
        }
    return {
        "schema": "legacy-qha-delivery/v1",
        "produced_by": "ht-vasp/scripts/build_legacy_qha_delivery.py (A line)",
        "read_this_first": "README.md",
        "layout_is_a_proposal_not_a_contract": "LAYOUT_PROPOSAL.md",
        "n_endpoints": len(endpoints),
        "endpoints": endpoints,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=None,
                    help="the endpoint root that holds legacy-as-received/ "
                         "(not needed with --index-only)")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--label", default=None)
    ap.add_argument("--external-structure-json", type=Path, default=None,
                    help="a pymatgen Structure JSON (e.g. an ht-vasp static.json) "
                         "whose relaxed cell is copied in as structure/from-ht-vasp")
    ap.add_argument("--external-structure-template", type=Path, default=None,
                    help="an ht-vasp template POSCAR to copy in as a labelled "
                         "reference (its lattice is a shared placeholder)")
    ap.add_argument("--json-out", type=Path, default=None,
                    help="write the endpoint summary here instead of stdout only")
    ap.add_argument("--content-dir", type=Path, default=None,
                    help="hand-authored documents to overlay (shared files at its "
                         "root, per-endpoint files under <content-dir>/<label>/)")
    ap.add_argument("--index-only", action="store_true",
                    help="skip building and just rewrite <outdir>/DELIVERY.json "
                         "from the endpoints already on disk")
    args = ap.parse_args()

    if args.index_only:
        idx = build_index(args.outdir)
        (args.outdir / "DELIVERY.json").write_text(
            json.dumps(idx, indent=2) + "\n", encoding="utf-8")
        print(f"DELIVERY.json: {idx['n_endpoints']} endpoints")
        for k, v in idx["endpoints"].items():
            print(f"  {k}: {v['audit_status']} "
                  f"({v['audit_n_failures']} failures, "
                  f"{v['audit_n_warnings']} warnings), "
                  f"applicability={v['applicability']}")
        return

    if args.root is None:
        ap.error("--root is required unless --index-only is given")
    label = args.label or args.root.name
    summary = build(args.root, args.outdir, label,
                    args.external_structure_json, args.external_structure_template)
    if args.content_dir is not None:
        generated = {(args.outdir / f["delivery_path"]).resolve()
                     for f in summary["files"]}
        generated |= {(args.outdir / label / "INVENTORY.json").resolve()}
        authored = overlay_content(args.content_dir, args.outdir / label, label,
                                   generated)
        summary["authored_documents"] = authored
    text = json.dumps(summary, indent=2) + "\n"
    (args.outdir / label / "INVENTORY.json").write_text(text, encoding="utf-8")
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(text, encoding="utf-8")
    bad = [f for f in summary["files"] if f["matches_mirror_manifest"] is False]
    print(f"{label}: {len(summary['files'])} files copied, "
          f"{len(bad)} failing the mirror-manifest hash check")
    print(f"  audit status: {summary['audit']['status']} "
          f"({summary['audit']['n_failures']} failures, "
          f"{summary['audit']['n_warnings']} warnings)")
    print(f"  raw VASP layer: {summary['raw_index']['n_files']} files, "
          f"{summary['raw_index']['total_bytes'] / 2**20:.2f} MiB (referenced, not copied)")
    for f in bad:
        print(f"  HASH MISMATCH {f['mirror_path']}")


if __name__ == "__main__":
    main()
