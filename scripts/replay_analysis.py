#!/usr/bin/env python3
"""Re-derive an endpoint's stability verdict and QHA analysis from its store.

Why this exists: the stability verdict is a pure function of the stored force
constants, the cell and the thresholds, but the verdict is written into
``runs/<run_id>/<endpoint>-dfpt-qha.json``, which the flow treats as an
immutable artifact.  Changing a threshold therefore used to mean re-running
seven DFPT volumes -- hours of VASP -- to recompute a decision that needs no
VASP at all.

``collect validated volumes`` in the store already holds every volume's force
constants, cell and supercell, so the three downstream stages can be replayed:

    collect validated volumes
      -> classify_dfpt_stability   (the threshold being changed lives here)
      -> dfpt_thermal_properties
      -> analyze_dfpt_qha

The replay writes a NEW analysis file beside the old one and records both
hashes in a sidecar; it never edits or deletes the original.  That keeps
"what the run produced" and "what the current thresholds say" separable, which
is the whole point of an immutable run directory.

    python scripts/replay_analysis.py FCC-Nb-Ta-4 --dry-run
    python scripts/replay_analysis.py FCC-Nb-Ta-4
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

POOL_RECIPES = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json")

#: Same values main_dfpt_qha.py uses; recorded in the sidecar so a replay is
#: reproducible even if the defaults change later.
SUPERCELL_MATRIX = ((2, 0, 0), (0, 2, 0), (0, 0, 2))
MESH = (12, 12, 12)
TEMPERATURES = list(range(0, 1001, 100))
PRESSURE_GPA = 0.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def latest_run(endpoint: str) -> Path:
    """The run whose store is worth replaying.

    An endpoint can have several run directories (a failed attempt followed by
    a successful one).  Sorting by name does not order them by time -- the run
    id is a random suffix -- so this prefers a run that actually produced
    volumes, and breaks ties by modification time.
    """
    runs = list((ROOT / "data" / "endpoints" / endpoint / "dfpt-qhaflow" / "runs")
                .glob("dfpt-qha-*"))
    if not runs:
        raise SystemExit(f"{endpoint}: no runs under data/endpoints")

    def score(run: Path):
        store = run / "store.json"
        if not store.is_file():
            return (0, 0.0)
        try:
            records = json.loads(store.read_text())
        except Exception:
            return (0, 0.0)
        has_volumes = any((r.get("name") or "") == "collect validated volumes"
                          for r in records)
        return (1 if has_volumes else 0, store.stat().st_mtime)

    return max(runs, key=score)


def store_records(run: Path) -> list[dict]:
    return json.loads((run / "store.json").read_text())


def recorrect_force_constants(collected: dict, *, threshold: float) -> tuple[dict, list[dict]]:
    """Re-apply the ASR correction to stored matrices under the current limit.

    The correction runs in the collector, so what ``collect validated volumes``
    stores is the matrix *as the collector left it*.  A volume whose residual
    exceeded the old limit is therefore stored uncorrected, and changing the
    limit has to be applied here -- the classification stage only reads the
    matrix it is handed.
    """
    import numpy as np
    from phonopy.harmonic.force_constants import (
        set_permutation_symmetry, set_translational_invariance,
    )

    from htvasp.dfpt_qha.stability import ASR_TOL_REL

    applied: list[dict] = []
    volumes = []
    for record in collected.get("volumes") or []:
        fc = np.array(record["force_constants"], dtype=float)
        scale = float(np.abs(fc).max()) or 1.0
        before = float(np.abs(fc.sum(axis=1)).max()) / scale
        changed = False
        if before > ASR_TOL_REL and before <= threshold:
            corrected = np.array(fc, dtype=float)
            set_translational_invariance(corrected)
            set_permutation_symmetry(corrected)
            after = float(np.abs(corrected.sum(axis=1)).max()) / scale
            if after <= ASR_TOL_REL:
                fc = corrected
                changed = True
                applied.append({"index": record["index"], "asr_rel_before": before,
                                "asr_rel_after": after})
        volumes.append({**record, "force_constants": fc.tolist()})
    return {**collected, "volumes": volumes}, applied


def load_recipe(endpoint: str):
    from htvasp.model.dfpt_qha_recipe import Recipe

    recipes = json.loads(POOL_RECIPES.read_text())["endpoints"]
    if endpoint not in recipes:
        raise SystemExit(f"{endpoint} is not in {POOL_RECIPES}")
    return Recipe.from_dict(recipes[endpoint]["recipe"])


def replay(endpoint: str, *, write: bool) -> int:
    from htvasp.dfpt_qha.flow_jobs import classify_dfpt_stability
    from htvasp.dfpt_qha.qha import analyze_dfpt_qha, dfpt_thermal_properties
    from htvasp.dfpt_qha.stability import ASR_CORRECTION_MAX_REL

    run = latest_run(endpoint)
    records = store_records(run)
    collected = next((r["output"] for r in records
                      if (r.get("name") or "") == "collect validated volumes"), None)
    if collected is None:
        raise SystemExit(f"{run}: no 'collect validated volumes' record")

    old_json = run / f"{endpoint}-dfpt-qha.json"
    old = json.loads(old_json.read_text()) if old_json.is_file() else {}

    recipe = load_recipe(endpoint)
    collected, re_corrected = recorrect_force_constants(
        collected, threshold=ASR_CORRECTION_MAX_REL)

    print(f"{endpoint}  ({run.name})")
    print(f"  ASR_CORRECTION_MAX_REL = {ASR_CORRECTION_MAX_REL}")
    print(f"  volumes in store: {len(collected.get('volumes') or [])}")
    if re_corrected:
        for item in re_corrected:
            print(f"  re-corrected volume {item['index']}: "
                  f"{item['asr_rel_before']:.5f} -> {item['asr_rel_after']:.2e}")
    else:
        print("  no volume needed re-correction under the current limit")
    print(f"  old verdict: {old.get('qha_applicability')}  ok={old.get('ok')}")

    # These three are jobflow @job functions.  Calling one returns a Response
    # whose ``output`` is an OutputReference, not the value; the value only
    # exists once jobflow has run it.  So call the underlying function with
    # .original, which is the plain Python callable.
    classification = classify_dfpt_stability.original(
        collected, recipe=recipe.model_dump(mode="json"),
        supercell_matrix=SUPERCELL_MATRIX, mesh=MESH,
    ).output
    new_applicability = classification.get("qha_applicability")
    statuses = [r.get("status") for r in classification.get("reports") or []]
    print(f"  new verdict: {new_applicability}")
    print(f"  per-volume: {statuses}")

    thermal = dfpt_thermal_properties.original(
        classification, supercell_matrix=SUPERCELL_MATRIX, mesh=MESH,
        temperatures=TEMPERATURES,
    )
    response = analyze_dfpt_qha.original(thermal, pressure_gpa=PRESSURE_GPA)
    analysis = dict(response.output)
    analysis["endpoint"] = endpoint
    analysis["replayed_from"] = {
        "run": run.name,
        "store_sha256": sha256(run / "store.json"),
        "original_analysis_sha256": sha256(old_json) if old_json.is_file() else None,
        "asr_correction_max_rel": ASR_CORRECTION_MAX_REL,
        "volumes_re_corrected": re_corrected,
        "supercell_matrix": [list(r) for r in SUPERCELL_MATRIX],
        "mesh": list(MESH),
        "temperatures": TEMPERATURES,
        "pressure_gpa": PRESSURE_GPA,
    }

    print(f"  replayed ok={analysis.get('ok')}  "
          f"V0={analysis.get('eos', {}).get('V0_angstrom3')}")
    if not write:
        print("  dry-run: nothing written")
        return 0

    out_json = run / f"{endpoint}-dfpt-qha.replayed.json"
    out_json.write_text(json.dumps(analysis, indent=2, sort_keys=True) + "\n")
    sidecar = run / f"{endpoint}-dfpt-qha.replay.json"
    sidecar.write_text(json.dumps({
        "schema": "htvasp-dfpt-qha-replay/v1",
        "endpoint": endpoint,
        "run": run.name,
        "original": {"path": old_json.name,
                     "qha_applicability": old.get("qha_applicability"),
                     "sha256": analysis["replayed_from"]["original_analysis_sha256"]},
        "replayed": {"path": out_json.name,
                     "qha_applicability": new_applicability,
                     "sha256": sha256(out_json)},
        "reason": "threshold change; no VASP re-run",
    }, indent=2, sort_keys=True) + "\n")
    print(f"  wrote {out_json.name} and {sidecar.name}")
    print("  the original analysis file is untouched")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("endpoints", nargs="+")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    failures = 0
    for endpoint in args.endpoints:
        try:
            replay(endpoint, write=not args.dry_run)
        except SystemExit as exc:
            print(f"{exc}", file=sys.stderr)
            failures += 1
        print()
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
