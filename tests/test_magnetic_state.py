"""Tests for reading the converged magnetic state out of an archived OUTCAR.

The point of the script is a distinction that a naive reader gets wrong: a net
moment of zero is not the same as a nonmagnetic cell.  These tests pin that
down, plus the OUTCAR parsing that feeds it.
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import extract_magnetic_state as M  # noqa: E402

BLOCK = """ magnetization (x)

# of ion       s       p       d       tot
------------------------------------------
{rows}
--------------------------------------------------
tot          {tot}
"""

def _outcar(tmp_path: Path, moments: list[float], extra_blocks: int = 0) -> Path:
    """An OUTCAR holding ``extra_blocks`` earlier blocks then ``moments``."""
    def block(values):
        rows = "\n".join(f"    {i + 1}      0.01   0.02   0.00  {v:8.3f}"
                         for i, v in enumerate(values))
        return BLOCK.format(rows=rows, tot=f"{sum(values):8.3f}")
    text = "".join(block([9.0] * len(moments)) for _ in range(extra_blocks))
    text += block(moments)
    path = tmp_path / "OUTCAR.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(text)
    return path

def test_last_block_wins(tmp_path):
    """The final ionic step is the answer; earlier ones are history."""
    path = _outcar(tmp_path, [1.1, -1.1], extra_blocks=2)
    assert M.per_site_moments(path) == pytest.approx([1.1, -1.1])

def test_zero_net_with_large_abs_is_antiferromagnetic():
    record = M.classify([1.134, -1.134])
    assert record["kind"] == "antiferromagnetic"
    assert record["net_moment_muB"] == pytest.approx(0.0, abs=1e-9)
    assert record["mean_abs_muB"] == pytest.approx(1.134)
    # IHJ still needs a magnitude, and the sign of the ordering lives on TC
    assert record["bmagn_ihj_muB"] == pytest.approx(1.134)
    assert record["tc_sign"] == -1

def test_zero_net_with_zero_abs_is_nonmagnetic():
    record = M.classify([0.0, -0.0])
    assert record["kind"] == "nonmagnetic"
    assert record["tc_sign"] == 0

def test_ferromagnetic_keeps_its_sign_convention():
    record = M.classify([2.5, -0.04])
    assert record["kind"] == "ferromagnetic"
    assert record["tc_sign"] == 1
    assert record["bmagn_ihj_muB"] == pytest.approx((2.5 + 0.04) / 2)

def test_initial_moments_are_read_from_the_incar(tmp_path):
    path = _outcar(tmp_path, [1.0, 2.0])
    with gzip.open(tmp_path / "INCAR.gz", "wt") as handle:
        handle.write("ALGO = Normal\nMAGMOM = 1*1.0 1*3.0\n")
    assert M.initial_moments(path) == "1*1.0 1*3.0"

def test_missing_magnetization_block_is_reported_not_guessed(tmp_path):
    path = tmp_path / "OUTCAR.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("no magnetization here\n")
    assert M.per_site_moments(path) is None

def test_a_collapsed_magnetic_run_is_flagged(tmp_path, monkeypatch):
    """A recipe that expects magnetism but gets 0 mu_B must not pass silently."""
    run = (tmp_path / "data" / "endpoints" / "BCC-X-X-2" / "dfpt-qhaflow"
           / "runs" / "r" / "1-EOS_equilibrium_relaxation")
    run.mkdir(parents=True)
    _outcar(run, [0.0, 0.0])
    with gzip.open(run / "INCAR.gz", "wt") as handle:
        handle.write("MAGMOM = 2*3.0\n")

    pool = tmp_path / "pool.json"
    pool.write_text(json.dumps({"endpoints": {"BCC-X-X-2": {
        "recipe": {"magnetic_policy": "confirmed_branch"}}}}))
    monkeypatch.setattr(M, "POOL", pool)
    monkeypatch.setattr(M, "ROOT", tmp_path)

    assert M.main(["-n", "BCC-X-X-2"]) == 0
    record = json.loads((tmp_path / "data" / "endpoints" / "BCC-X-X-2"
                         / "magnetic-state.json").read_text())
    assert record["kind"] == "nonmagnetic"
    assert "expects a magnetic branch" in record["warning"]
