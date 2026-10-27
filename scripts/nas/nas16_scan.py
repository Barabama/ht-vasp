#!/usr/bin/env python
"""NAS-side fail-closed reliability scan over the legacy DFPT-QHA endpoints.

Runs ON the NAS (stdlib only, Python 3.14). Never writes outside the audit tree
it is told to use. Emits one JSON object per endpoint plus a summary JSON.

What it checks, per endpoint, from the files that are actually there:

  R3 layer      ISIF/IBRION/NSW/EDIFF/ISMEAR, k-mesh, species order vs POTCAR,
                the "reached required accuracy" marker, the ICHARG=1 charge
                density story (incomplete / foreign / dimension mismatch /
                atomic superposition), restart-log count
  DFPT layer    IBRION=8 present, number of symmetry-reduced perturbations,
                whether every one was solved, the ground-state SCF loop shape
                (iterations vs NELM, tail energy span), the DFPT ground-state
                energy against the same volume's R3f static
  FORCE_CONSTANTS   full (n,n,3,3) parse: max |Phi|, acoustic-sum-rule residual,
                permutation-symmetry deviation
  R3f statics   convergence markers and energy
  thermal yaml  num_modes vs num_integrated_modes (phonopy drops imaginary modes
                silently), ZPE, high-T entropy, temperature grid
  v-e.dat       row count, spacing, whether row i equals volume i's R3f energy
  curves        Gibbs temperature grid, largest adjacent jump
  cross-consistency  same-geometry energy differences between R3, R3f and the
                DFPT ground state; whether the curve's grid matches the cells

Usage (on the NAS):
    python scan.py --root "E:\\nas-shared\\...\\16symbols-cxy-hmd" \
                   --out  "E:\\nas-shared\\...\\_reliability_audit\\data" \
                   [--only BCC-Co-Ti-2] [--limit 20]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import time

# ------------------------------------------------------------------ helpers --

ENERGY_RE = re.compile(
    r"energy\s+without entropy=\s*(-?[\d.]+)\s+energy\(sigma->0\)\s*=\s*(-?[\d.]+)")
DOF_RE = re.compile(r"Degree of freedom:\s*(\d+)\s*/\s*(\d+)")
ITER_RE = re.compile(r"^\s*(?:DAV|RMM):\s*\d+\s+(-?[\d.]+E?[-+]?\d*)\s+(-?[\d.]+E?[-+]?\d*)",
                     re.M)
VOL_RE = re.compile(r"^QHA-(\d+)$")

COMPLETION = "General timing and accounting informations for this job"
SCF_OK = "aborting loop because EDIFF is reached"
SCF_BAD = "aborting loop EDIFF was not reached (unconverged)"
ATOMIC_GUESS = "charge density of overlapping atoms calculated"


def read(path, limit=None):
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            return f.read() if limit is None else f.read(limit)
    except OSError:
        return None


def num(path):
    """Read a whitespace-separated numeric table into a list of float rows."""
    txt = read(path)
    if txt is None:
        return None
    rows = []
    for line in txt.splitlines():
        s = line.split()
        if not s or s[0].startswith("#"):
            continue
        try:
            rows.append([float(x) for x in s])
        except ValueError:
            continue
    return rows


def parse_fc(path):
    """Full (n, n, 3, 3) FORCE_CONSTANTS. Returns None if unreadable."""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            head = [int(x) for x in f.readline().split()]
            n, m = (head[0], head[0]) if len(head) == 1 else (head[0], head[1])
            fc = [[[[0.0] * 3 for _ in range(3)] for _ in range(m)] for _ in range(n)]
            for i in range(n):
                for j in range(m):
                    f.readline()
                    for k in range(3):
                        v = [float(x) for x in f.readline().split()]
                        fc[i][j][k] = v
    except (OSError, ValueError, IndexError):
        return None
    return fc


def fc_stats(fc):
    n = len(fc)
    m = len(fc[0])
    maxabs = 0.0
    for i in range(n):
        for j in range(m):
            for k in range(3):
                for l in range(3):
                    a = abs(fc[i][j][k][l])
                    if a > maxabs:
                        maxabs = a
    # ASR: sum_j Phi[i][j] should vanish for every i, k, l
    asr = 0.0
    for i in range(n):
        for k in range(3):
            for l in range(3):
                s = sum(fc[i][j][k][l] for j in range(m))
                asr = max(asr, abs(s))
    # permutation symmetry Phi[i][j][k][l] == Phi[j][i][l][k]
    perm = 0.0
    if n == m:
        for i in range(n):
            for j in range(n):
                for k in range(3):
                    for l in range(3):
                        perm = max(perm, abs(fc[i][j][k][l] - fc[j][i][l][k]))
    return {"shape": [n, m, 3, 3], "max_abs": maxabs, "asr_abs": asr,
            "asr_relative": (asr / maxabs if maxabs else None),
            "permutation_max_dev": perm}


def scf_shape(text, upto_marker):
    """Iteration count and tail wander of the SCF block before a marker."""
    if text is None:
        return {}
    cut = text.find(upto_marker) if upto_marker else -1
    block = text[:cut] if cut > 0 else text
    it = ITER_RE.findall(block)
    out = {"iterations": len(it)}
    if it:
        tail = [float(e) for _, e in it[-10:]]
        out["last_dE_eV"] = float(it[-1][1])
        out["tail_span_meV"] = (max(tail) - min(tail)) * 1000.0
    return out


def energy_blocks(text):
    if text is None:
        return []
    return ENERGY_RE.findall(text)


def yaml_scalars(path):
    """Minimal reader for the flat phonopy thermal_properties yaml header."""
    txt = read(path)
    if txt is None:
        return None
    out = {}
    for line in txt.splitlines():
        if line.startswith("- temperature:"):
            break
        m = re.match(r"^([a-zA-Z_0-9]+):\s*(\S+)\s*$", line)
        if m:
            k, v = m.group(1), m.group(2)
            try:
                out[k] = float(v)
            except ValueError:
                out[k] = v
    return out


def yaml_temperature_grid(path):
    txt = read(path)
    if txt is None:
        return None
    ts = [float(m.group(1)) for m in
          re.finditer(r"^-\s*temperature:\s*(-?[\d.]+)", txt, re.M)]
    if not ts:
        return None
    return {"n": len(ts), "first": ts[0], "last": ts[-1],
            "step": (ts[1] - ts[0]) if len(ts) > 1 else None}


def poscar_lattice(path):
    txt = read(path)
    if txt is None:
        return None
    lines = [l for l in txt.splitlines()]
    if len(lines) < 5:
        return None
    try:
        scale = float(lines[1].split()[0])
        lat = [[float(x) * scale for x in lines[2 + i].split()[:3]] for i in range(3)]
    except (ValueError, IndexError):
        return None
    det = (lat[0][0] * (lat[1][1] * lat[2][2] - lat[1][2] * lat[2][1])
           - lat[0][1] * (lat[1][0] * lat[2][2] - lat[1][2] * lat[2][0])
           + lat[0][2] * (lat[1][0] * lat[2][1] - lat[1][1] * lat[2][0]))
    return {"a": math.sqrt(sum(x * x for x in lat[0])), "volume": abs(det)}


def incar_tags(path):
    txt = read(path)
    if txt is None:
        return None
    out = {}
    for line in txt.splitlines():
        line = line.split("!")[0].split("#")[0]
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip().upper()] = v.strip().split()[0] if v.strip() else ""
    return out


def potcar_species(path):
    txt = read(path)
    if txt is None:
        return None
    return [m.group(1) for m in re.finditer(r"TITEL\s*=\s*PAW\w*\s+(\w+)", txt)]


def poscar_species(path):
    """The species line of a POSCAR (the real one, not the free-text title)."""
    txt = read(path)
    if txt is None:
        return None
    lines = txt.splitlines()
    for i in range(5, min(9, len(lines))):
        s = lines[i].split()
        if len(s) >= 2 and all(re.fullmatch(r"[A-Za-z]+", x) for x in s):
            return s
    return None


def poscar_title(path):
    txt = read(path)
    if txt is None:
        return None
    lines = txt.splitlines()
    return lines[0].strip() if lines else None


def poscar_counts(path):
    """The per-species atom counts, and their sum."""
    txt = read(path)
    if txt is None:
        return None, None
    lines = txt.splitlines()
    start = None
    for i in range(5, min(9, len(lines))):
        s = lines[i].split()
        if len(s) >= 2 and all(re.fullmatch(r"[A-Za-z]+", x) for x in s):
            start = i + 1
            break
    if start is None:
        for i in range(5, min(9, len(lines))):
            s = lines[i].split()
            if len(s) >= 1 and all(re.fullmatch(r"\d+", x) for x in s):
                start = i
                break
    if start is None or start >= len(lines):
        return None, None
    s = lines[start].split()
    if not all(re.fullmatch(r"\d+", x) for x in s):
        return None, None
    c = [int(x) for x in s]
    return c, sum(c)


def parse_v_e(path):
    rows = num(path)
    if not rows:
        return None
    return {"n": len(rows),
            "volumes": [r[0] for r in rows],
            "energies": [r[1] for r in rows],
            "step": (rows[1][0] - rows[0][0]) if len(rows) > 1 else None}


def gibbs_jump(path):
    rows = num(path)
    if not rows or len(rows) < 2:
        return None
    ts = [r[0] for r in rows]
    gs = [r[1] for r in rows]
    d = [abs(gs[i + 1] - gs[i]) for i in range(len(gs) - 1)]
    i = max(range(len(d)), key=lambda k: d[k])
    return {"n": len(rows), "first_T": ts[0], "last_T": ts[-1],
            "step_T": ts[1] - ts[0],
            "G_first": gs[0], "G_last": gs[-1],
            "max_jump_eV": d[i], "max_jump_at_K": [ts[i], ts[i + 1]],
            "n_jumps_over_0.2": sum(1 for x in d if x > 0.2)}


# ------------------------------------------------------------------- scan ----

def scan_endpoint(root, ep):
    d = os.path.join(root, ep)
    r = {"endpoint": ep, "findings": [], "warnings": [], "info": {}}

    def fail(msg):
        r["findings"].append(msg)

    def warn(msg):
        r["warnings"].append(msg)

    def listdir(p):
        try:
            return sorted(os.listdir(p))
        except OSError:
            return None

    top = listdir(d)
    if top is None:
        fail("endpoint directory unreadable")
        return r

    # An endpoint with no VASP layer cannot be verified at all. Report that as
    # its own state rather than as a pile of missing-file failures.
    has_dfpt_dir = os.path.isdir(os.path.join(d, "dfpt"))
    has_r3_dir = (os.path.isdir(os.path.join(d, "R3"))
                  or os.path.isdir(os.path.join(d, "R31")))
    r["info"]["vasp_layer_present"] = bool(has_dfpt_dir or has_r3_dir)

    # ---- recipe at the top level ---------------------------------------
    inc_r31 = incar_tags(os.path.join(d, "INCAR-R31"))
    inc_r3f = incar_tags(os.path.join(d, "INCAR-R3f"))
    inc_ph = incar_tags(os.path.join(d, "INCAR-Phonopy"))
    r["info"]["recipe"] = {"R31": inc_r31, "R3f": inc_r3f, "Phonopy": inc_ph}
    if inc_r31 is None and r["info"]["vasp_layer_present"]:
        fail("no INCAR-R31 at the endpoint root")
    elif inc_r31 is not None:
        if inc_r31.get("ISIF") != "3":
            warn(f"INCAR-R31 ISIF={inc_r31.get('ISIF')} (expected 3)")
        if inc_r31.get("IBRION") != "2":
            warn(f"INCAR-R31 IBRION={inc_r31.get('IBRION')} (expected 2 for the ionic relax)")
    if inc_ph is not None and inc_ph.get("IBRION") != "8":
        fail(f"INCAR-Phonopy IBRION={inc_ph.get('IBRION')} (expected 8 for DFPT)")

    pot = potcar_species(os.path.join(d, "POTCAR"))
    pos = poscar_species(os.path.join(d, "POSCAR"))
    title = poscar_title(os.path.join(d, "POSCAR"))
    r["info"]["potcar_species"] = pot
    r["info"]["poscar_species"] = pos
    r["info"]["poscar_title"] = title
    if pot and pos and len(pot) > len(pos):
        # VASP reads only as many POTCAR blocks as POSCAR declares species, and
        # it uses POTCAR's order, silently overwriting POSCAR's names. So what
        # matters is that the first len(pos) blocks ARE the intended pair.
        warn(f"POTCAR has {len(pot)} blocks ({pot}) but POSCAR declares "
             f"{len(pos)} species ({pos}); VASP reads the first {len(pos)}")
    if pot and pos:
        head = pot[:len(pos)]
        if [x.upper() for x in head] != [x.upper() for x in pos]:
            fail(f"POTCAR's first {len(pos)} blocks are {head} but POSCAR declares "
                 f"{pos}: VASP silently uses the POTCAR order, so the elements are "
                 f"attached to the wrong sites")
    if title and pos and title.split() != pos:
        r["info"]["poscar_title_is_not_the_species_line"] = True

    # ---- R3 layer ------------------------------------------------------
    r3 = os.path.join(d, "R3")
    if not os.path.isdir(r3):
        r3alt = os.path.join(d, "R31")
        if os.path.isdir(r3alt):
            r["info"]["r3_dir"] = "R31"
            r3 = r3alt
        else:
            r["info"]["r3_dir"] = None
    else:
        r["info"]["r3_dir"] = "R3"

    if os.path.isdir(r3):
        r3files = listdir(r3) or []
        r["info"]["r3_files"] = r3files
        out = read(os.path.join(r3, "OUTCAR"))
        if out is None:
            fail("R3/OUTCAR missing")
        else:
            r["info"]["r3_reached_required_accuracy"] = out.count("reached required accuracy")
            r["info"]["r3_scf_ok"] = out.count(SCF_OK)
            r["info"]["r3_scf_bad"] = out.count(SCF_BAD)
            r["info"]["r3_completion"] = COMPLETION in out
            r["info"]["r3_ionic_steps"] = out.count("free  energy   TOTEN")
            r["info"]["r3_atomic_guess"] = ATOMIC_GUESS in out
            eb = energy_blocks(out)
            r["info"]["r3_last_energy_sigma0"] = float(eb[-1][1]) if eb else None
            if r["info"]["r3_reached_required_accuracy"] == 0:
                fail("R3 relaxation never printed 'reached required accuracy'")
            if not r["info"]["r3_completion"]:
                fail("R3/OUTCAR has no job-completion marker")
            if r["info"]["r3_scf_bad"] > 0:
                fail(f"R3/OUTCAR reports {r['info']['r3_scf_bad']} unconverged SCF loop(s)")
        # restart logs
        logs = [f for f in r3files if re.fullmatch(r"R\d+", f)]
        r["info"]["r3_restart_logs"] = logs
        stories = {}
        for lg in logs:
            t = read(os.path.join(r3, lg)) or ""
            st = []
            if "chargedensity file is incomplete" in t:
                st.append("incomplete_chargedensity")
            m = re.search(r"charge-density read from file:\s*(\S+)", t)
            if m:
                st.append("read_foreign_chargedensity:" + m.group(1))
            if "dimensions on CHGCAR file are different" in t:
                st.append("chgcar_dimension_mismatch")
            if ATOMIC_GUESS in t:
                st.append("atomic_superposition")
            if not st:
                st.append("clean")
            stories[lg] = st
        r["info"]["r3_restart_stories"] = stories
        bad_stories = {k: v for k, v in stories.items() if v != ["clean"]}
        if bad_stories:
            warn(f"R3 relaxation had {len(bad_stories)} non-clean restart(s): {bad_stories}")

        # k-mesh
        kp = read(os.path.join(r3, "KPOINTS")) or ""
        r["info"]["r3_kpoints"] = [l.strip() for l in kp.splitlines() if l.strip()][:6]

    # ---- curve location -------------------------------------------------
    # Two layouts exist: flat (curves at the endpoint root, e.g. every FCC and
    # BCC-Nb-Ti-2) and nested (curves inside QHA-<name>/, e.g. BCC-Al-Co-2).
    curve_dir = None
    for name in listdir(d) or []:
        p = os.path.join(d, name)
        if os.path.isdir(p) and name not in ("dfpt", "R3", "R31") and name != "_reliability_audit":
            if os.path.isfile(os.path.join(p, "gibbs-temperature.dat")):
                curve_dir = name
                break
    r["info"]["curve_dir"] = curve_dir
    r["info"]["curve_layout"] = "nested" if curve_dir else (
        "flat" if os.path.isfile(os.path.join(d, "gibbs-temperature.dat")) else None)
    curve_root = os.path.join(d, curve_dir) if curve_dir else d
    ve = parse_v_e(os.path.join(d, "v-e.dat"))
    if ve is None:
        ve = parse_v_e(os.path.join(curve_root, "v-e.dat"))
    r["info"]["v_e"] = ve
    if ve is None:
        fail("no readable v-e.dat")

    # ---- is the E(V) table actually a curve? ----------------------------
    # A QHA volume scan spans a few percent of the volume, so its energies must
    # span a small range and be smooth. An energy span of electronvolts, or a
    # second difference comparable to the span, means the static calculations
    # did not converge to one curve at all.
    if ve and ve["n"] >= 4:
        E = ve["energies"]
        span = max(E) - min(E)
        rough = max(abs(E[i + 1] - 2 * E[i] + E[i - 1])
                    for i in range(1, len(E) - 1))
        i0 = min(range(len(E)), key=lambda k: E[k])
        r["info"]["v_e_shape"] = {
            "span_eV": span,
            "max_second_difference_eV": rough,
            "roughness_ratio": rough / span if span else None,
            "min_at_index": i0,
            "min_is_interior": 0 < i0 < len(E) - 1,
        }
        if span > 1.0:
            fail(f"v-e.dat spans {span:.3f} eV across the volume scan: these are "
                 f"not converged points of one curve")
        elif rough > 0.05 and rough / span > 0.5:
            fail(f"v-e.dat is not smooth: a second difference of {rough:.4f} eV "
                 f"against a total span of {span:.4f} eV")
        elif span > 0.5:
            warn(f"v-e.dat spans {span:.3f} eV, much wider than a QHA scan should")
    gj = gibbs_jump(os.path.join(curve_root, "gibbs-temperature.dat"))
    r["info"]["gibbs"] = gj
    if gj and gj["max_jump_eV"] > 0.2:
        warn(f"reported Gibbs has a {gj['max_jump_eV']:.4f} eV jump at "
             f"{gj['max_jump_at_K']} and {gj['n_jumps_over_0.2']} jumps over 0.2 eV")

    # ---- every thermal_properties yaml we can see, wherever it lives ---
    # Flat endpoints keep them at the root; the nested ones keep them in both
    # dfpt/QHA-NN/ and the curve directory. Collect by volume index.
    tp_files = {}
    for f in (top or []):
        m = re.match(r"^thermal_properties-(\d+)\.yaml$", f)
        if m:
            tp_files.setdefault(int(m.group(1)), []).append(os.path.join(d, f))
    if curve_dir:
        for f in (listdir(curve_root) or []):
            m = re.match(r"^thermal_properties-(\d+)\.yaml$", f)
            if m:
                tp_files.setdefault(int(m.group(1)), []).append(
                    os.path.join(curve_root, f))
    tp_info = {}
    for idx in sorted(tp_files):
        s = yaml_scalars(tp_files[idx][0]) or {}
        g = yaml_temperature_grid(tp_files[idx][0])
        tp_info[idx] = {"file": os.path.relpath(tp_files[idx][0], d),
                        "num_modes": s.get("num_modes"),
                        "num_integrated_modes": s.get("num_integrated_modes"),
                        "natom": s.get("natom"),
                        "zpe_kJ_per_mol": s.get("zero_point_energy"),
                        "high_T_entropy": s.get("high_T_entropy"),
                        "grid": g}
    if tp_info:
        r["info"]["thermal_properties"] = tp_info
        dropped = [k for k, v in tp_info.items()
                   if v["num_modes"] and v["num_integrated_modes"]
                   and v["num_modes"] > v["num_integrated_modes"]]
        if dropped:
            warn(f"{len(dropped)}/{len(tp_info)} thermal_properties yaml drop "
                 f"imaginary modes (phonon frequencies < 0 are silently excluded)")

    # ---- the one independent consistency check that works everywhere ----
    # For a QHA, G(T=0) = min_V E(V) + ZPE(at the equilibrium volume). If the
    # curve, the E(V) table and the phonon layer came from one calculation,
    # these three numbers must agree. Nothing here needs the VASP layer.
    if ve and tp_info and gj:
        E = ve["energies"]
        i0 = min(range(len(E)), key=lambda k: E[k])
        natom = None
        for k in sorted(tp_info):
            if tp_info[k].get("natom"):
                natom = tp_info[k]["natom"]
                break
        zpe = None
        if i0 + 1 in tp_info and tp_info[i0 + 1].get("zpe_kJ_per_mol") is not None:
            zpe = tp_info[i0 + 1]["zpe_kJ_per_mol"]
        elif tp_info:
            nearest = min(tp_info, key=lambda k: abs(k - (i0 + 1)))
            zpe = tp_info[nearest].get("zpe_kJ_per_mol")
        if zpe is not None:
            KJ_PER_EV = 96.48533212331002
            pred = E[i0] + zpe / KJ_PER_EV
            dev_meV = (gj["G_first"] - pred) * 1000.0
            r["info"]["g0_vs_minE_plus_zpe"] = {
                "G_at_0K_eV": gj["G_first"],
                "min_E_eV": E[i0],
                "min_E_at_volume": ve["volumes"][i0],
                "min_E_is_interior": 0 < i0 < len(E) - 1,
                "zpe_kJ_per_mol": zpe,
                "zpe_volume_index": (i0 + 1),
                "predicted_G0_eV": pred,
                "deviation_meV": dev_meV,
            }
            if abs(dev_meV) > 20.0:
                warn(f"G(0) - [min E(V) + ZPE] = {dev_meV:.1f} meV: the curve, "
                     f"the E(V) table and the phonon layer do not close")
        if not (0 < i0 < len(E) - 1):
            warn("the E(V) minimum is at an endpoint of the volume grid, so the "
                 "grid does not bracket the equilibrium")
        if ve["step"] and len(ve["volumes"]) > 2:
            sp = [ve["volumes"][i + 1] - ve["volumes"][i]
                  for i in range(len(ve["volumes"]) - 1)]
            if (max(sp) - min(sp)) > 0.02:
                warn(f"the v-e.dat volume grid is not uniform (spacings "
                     f"{min(sp):.3f}..{max(sp):.3f} A^3)")

    # ---- top-level thermal yaml (flat layout) --------------------------
    tp_top = sorted(f for f in (top or []) if f.startswith("thermal_properties-"))
    if tp_top and not tp_info:
        rows = []
        for f in tp_top:
            s = yaml_scalars(os.path.join(d, f)) or {}
            g = yaml_temperature_grid(os.path.join(d, f))
            rows.append({"file": f, "num_modes": s.get("num_modes"),
                         "num_integrated_modes": s.get("num_integrated_modes"),
                         "zpe": s.get("zero_point_energy"),
                         "high_T_entropy": s.get("high_T_entropy"),
                         "grid": g})
        r["info"]["top_thermal_properties"] = rows

    # ---- DFPT layer ----------------------------------------------------
    dfpt = os.path.join(d, "dfpt")
    if not os.path.isdir(dfpt):
        r["info"]["has_dfpt"] = False
        r["info"]["status"] = "unverifiable"
        r["info"]["verdict"] = ("no VASP layer (no dfpt/ and no R3/ or R31/): the "
                                "curves and thermal properties can be read, but how "
                                "they were produced cannot be checked here")
        return r
    r["info"]["has_dfpt"] = True

    vols = sorted((v for v in (listdir(dfpt) or []) if VOL_RE.match(v)),
                  key=lambda v: int(VOL_RE.match(v).group(1)))
    r["info"]["dfpt_volumes"] = vols
    if len(vols) != 10:
        warn(f"dfpt has {len(vols)} volume directories (the two known endpoints have 10)")

    vol_rows = []
    for v in vols:
        vd = os.path.join(dfpt, v)
        tag = VOL_RE.match(v).group(1)
        row = {"tag": tag}
        # the per-volume cell
        cell = poscar_lattice(os.path.join(vd, "POSCAR-unitcell"))
        row["cell"] = cell
        out = read(os.path.join(vd, "OUTCAR"))
        log = read(os.path.join(vd, "log-DFPT"))
        inc = incar_tags(os.path.join(vd, "INCAR"))
        row["incar"] = inc
        if out is None:
            row["error"] = "no OUTCAR"
            vol_rows.append(row)
            continue
        dofs = DOF_RE.findall(out)
        row["dfpt_perturbations"] = [[int(a), int(b)] for a, b in dofs]
        row["dfpt_all_solved"] = bool(dofs) and dofs[-1][0] == dofs[-1][1]
        row["completion"] = COMPLETION in out
        row["scf_ok"] = out.count(SCF_OK)
        row["scf_bad"] = out.count(SCF_BAD)
        row["atomic_guess"] = ATOMIC_GUESS in out
        row["linear_response_finished"] = (
            ("Linear response finished" in out)
            or (log is not None and "Linear response finished" in log))
        eb = energy_blocks(out)
        row["n_energy_blocks"] = len(eb)
        if eb:
            row["dfpt_gs_sigma0"] = float(eb[0][1])
        gs = scf_shape(log if log is not None else out,
                       "Linear response reoptimize")
        row["dfpt_gs_scf"] = gs
        nelm = None
        m = re.search(r"NELM\s*=\s*(\d+)", out)
        if m:
            nelm = int(m.group(1))
        row["nelm"] = nelm
        row["dfpt_gs_hit_nelm"] = bool(nelm and gs.get("iterations", 0) >= nelm)
        if not row["completion"]:
            fail(f"{v}: DFPT OUTCAR has no completion marker")
        if row["scf_ok"] == 0:
            fail(f"{v}: DFPT OUTCAR never reached EDIFF")
        if not row["dfpt_all_solved"]:
            fail(f"{v}: DFPT did not solve every perturbation ({row['dfpt_perturbations']})")
        if not row["linear_response_finished"]:
            fail(f"{v}: DFPT has no 'Linear response finished' marker")
        if row["atomic_guess"]:
            fail(f"{v}: DFPT fell back to the atomic superposition (ICHARG=1, no usable CHGCAR)")

        # R3f static
        r3f = os.path.join(vd, "R3f")
        r3fo = read(os.path.join(r3f, "OUTCAR"))
        if r3fo is None:
            fail(f"{v}: no R3f/OUTCAR (the static that produced this v-e.dat row)")
        else:
            r3fe = energy_blocks(r3fo)
            row["r3f_energy_sigma0"] = float(r3fe[0][1]) if r3fe else None
            row["r3f_scf_ok"] = r3fo.count(SCF_OK)
            row["r3f_scf_bad"] = r3fo.count(SCF_BAD)
            row["r3f_completion"] = COMPLETION in r3fo
            row["r3f_atomic_guess"] = ATOMIC_GUESS in r3fo
            r3flog = read(os.path.join(r3f, "R3f"))
            row["r3f_scf"] = scf_shape(r3flog if r3flog is not None else r3fo, None)
            if row["r3f_scf_ok"] == 0:
                fail(f"{v}: R3f static never reached EDIFF")
            if row["r3f_atomic_guess"]:
                fail(f"{v}: R3f static fell back to the atomic superposition")
            # DFPT ground state vs R3f static, per atom
            if row.get("dfpt_gs_sigma0") is not None and row["r3f_energy_sigma0"] is not None:
                _, n_at = poscar_counts(os.path.join(vd, "POSCAR"))
                _, n_r3f = poscar_counts(os.path.join(r3f, "POSCAR"))
                if n_r3f is None:
                    n_r3f = poscar_counts(os.path.join(vd, "POSCAR-unitcell"))[1]
                row["n_supercell_atoms"] = n_at
                row["r3f_n_atoms"] = n_r3f
                if n_at and n_r3f:
                    dmev = (row["dfpt_gs_sigma0"] / n_at
                            - row["r3f_energy_sigma0"] / n_r3f) * 1000.0
                    row["dfpt_vs_r3f_meV_per_atom"] = dmev
                    if abs(dmev) > 100.0:
                        fail(f"{v}: DFPT ground state is {dmev:.1f} meV/atom away from "
                             f"the R3f static it should share a geometry with")
        if row["dfpt_gs_hit_nelm"]:
            span = row["dfpt_gs_scf"].get("tail_span_meV")
            fail(f"{v}: DFPT ground-state SCF ran to NELM={nelm} without EDIFF"
                 + (f" (tail still moves {span:.1f} meV)" if span else ""))
        elif row["dfpt_gs_scf"].get("last_dE_eV") is not None:
            # Did the loop actually satisfy EDIFF, or just stop early?
            ed = None
            m2 = re.search(r"EDIFF\s*=\s*([\d.eE+-]+)", out)
            if m2:
                ed = float(m2.group(1))
            row["dfpt_ediif_eV"] = ed
            if ed is not None and abs(row["dfpt_gs_scf"]["last_dE_eV"]) > ed * 10:
                warn(f"{v}: DFPT ground-state SCF stopped at iteration "
                     f"{row['dfpt_gs_scf']['iterations']} with |dE|="
                     f"{abs(row['dfpt_gs_scf']['last_dE_eV']):.2e} eV, more than "
                     f"10x EDIFF={ed}")

        # FORCE_CONSTANTS
        fcs = parse_fc(os.path.join(vd, "FORCE_CONSTANTS"))
        if fcs is None:
            fail(f"{v}: no readable FORCE_CONSTANTS")
        else:
            st = fc_stats(fcs)
            row["force_constants"] = st
            if st["max_abs"] > 50.0:
                fail(f"{v}: FORCE_CONSTANTS max |Phi| = {st['max_abs']:.2f} eV/A^2 "
                     f"is not physical")
            if st["asr_relative"] and st["asr_relative"] > 0.10:
                fail(f"{v}: FORCE_CONSTANTS violates the acoustic sum rule by "
                     f"{st['asr_relative'] * 100:.1f}%")

        # thermal properties
        tp = [f for f in (listdir(vd) or []) if f.startswith("thermal_properties-")]
        if tp:
            s = yaml_scalars(os.path.join(vd, tp[0])) or {}
            g = yaml_temperature_grid(os.path.join(vd, tp[0]))
            row["thermal_properties"] = {
                "file": tp[0], "num_modes": s.get("num_modes"),
                "num_integrated_modes": s.get("num_integrated_modes"),
                "zpe": s.get("zero_point_energy"),
                "high_T_entropy": s.get("high_T_entropy"), "grid": g}
            if (s.get("num_modes") and s.get("num_integrated_modes")
                    and s["num_modes"] > s["num_integrated_modes"]):
                warn(f"{v}: thermal_properties drops "
                     f"{s['num_modes'] - s['num_integrated_modes']} modes "
                     f"(imaginary; phonopy does this silently)")

        # v-e.dat row vs this volume's R3f energy
        if ve and row.get("r3f_energy_sigma0") is not None:
            idx = int(tag) - 1
            if 0 <= idx < ve["n"]:
                dev = (ve["energies"][idx] - row["r3f_energy_sigma0"]) * 1000.0
                row["v_e_row_vs_r3f_meV"] = dev
                if abs(dev) > 1.0:
                    fail(f"{v}: v-e.dat row {idx} differs from its R3f energy by "
                         f"{dev:.3f} meV")
        vol_rows.append(row)

    r["info"]["volumes"] = vol_rows

    # ---- cross-checks --------------------------------------------------
    if ve and vol_rows and any(x.get("cell") for x in vol_rows):
        cells = [x.get("cell") for x in vol_rows]
        vs = [c["volume"] for c in cells if c]
        if len(vs) == len(ve["volumes"]):
            dev = max(abs(vs[i] - ve["volumes"][i]) for i in range(len(vs)))
            r["info"]["cell_vs_v_e_max_dev_ang3"] = dev
            if dev > 0.01:
                fail(f"per-volume cells disagree with v-e.dat volumes by up to "
                     f"{dev:.4f} A^3")

    if r["findings"]:
        r["info"]["status"] = "fail"
    elif r["info"].get("vasp_layer_present"):
        r["info"]["status"] = "pass"
    else:
        r["info"]["status"] = "unverifiable"
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", action="append", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--skip-no-dfpt", action="store_true",
                    help="only scan endpoints that carry a dfpt/ tree")
    args = ap.parse_args()

    root, out = args.root, args.out
    os.makedirs(out, exist_ok=True)

    if args.only:
        eps = args.only
    else:
        eps = sorted(x for x in os.listdir(root)
                     if os.path.isdir(os.path.join(root, x)) and not x.startswith("_"))
    if args.skip_no_dfpt:
        eps = [e for e in eps if os.path.isdir(os.path.join(root, e, "dfpt"))]
    if args.limit:
        eps = eps[:args.limit]

    results, t0 = [], time.time()
    for i, ep in enumerate(eps, 1):
        try:
            results.append(scan_endpoint(root, ep))
        except Exception as ex:                      # noqa: BLE001
            results.append({"endpoint": ep, "findings": [f"scan crashed: {ex!r}"],
                            "warnings": [], "info": {}})
        if i % 10 == 0 or i == len(eps):
            print(f"  {i}/{len(eps)}  ({time.time() - t0:.0f}s)", flush=True)

    summary = {"schema": "nas16-reliability-scan/v1", "root": root,
               "scanned_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "n_scanned": len(results),
               "n_with_findings": sum(1 for r in results if r["findings"]),
               "status_counts": {
                   s: sum(1 for r in results
                          if r["info"].get("status") == s)
                   for s in ("pass", "fail", "unverifiable")},
               "endpoints": results}
    p = os.path.join(out, "scan.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1)
    print(f"-> {p}")
    print(f"scanned {len(results)}, with findings "
          f"{summary['n_with_findings']}, in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
