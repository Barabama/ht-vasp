"""Original-versus-actual VASP input audit.

custodian may correct a run; what it rewrote has to be visible afterwards.  A
DFPT force-constant matrix is indexed by the POSCAR atom order, so a POSCAR
correction is the one that would silently mis-index everything downstream.
"""

from __future__ import annotations

import gzip

from htvasp.dfpt_qha.input_audit import audit_vasp_inputs

POSCAR = """B2 AlTa
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Al Ta
1 1
Direct
0.0 0.0 0.0
0.5 0.5 0.5
"""

KPOINTS = """Automatic
0
Gamma
8 8 8
0.0 0.0 0.0
"""


def _write(directory, name, text, gz=False):
    path = directory / (f"{name}.gz" if gz else name)
    if gz:
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return path


def _clean(directory, gz=False):
    _write(directory, "INCAR", "ALGO = Normal\nISMEAR = 1\nSIGMA = 0.2\n", gz)
    _write(directory, "KPOINTS", KPOINTS, gz)
    _write(directory, "POSCAR", POSCAR, gz)
    _write(directory, "POTCAR", "PAW_PBE Al\nPAW_PBE Ta\n", gz)
    for name in ("INCAR", "KPOINTS", "POSCAR", "POTCAR"):
        _write(directory, f"{name}.orig", _read(directory, name, gz), gz)


def _read(directory, name, gz=False):
    path = directory / (f"{name}.gz" if gz else name)
    if gz:
        with gzip.open(path, "rt") as handle:
            return handle.read()
    return path.read_text()


def test_clean_inputs_are_not_flagged(tmp_path):
    _clean(tmp_path)
    audit = audit_vasp_inputs(tmp_path)
    assert audit["modified_files"] == []
    assert audit["warnings"] == []
    assert audit["original_input_available"] is True
    assert audit["missing_original_files"] == []


def test_an_incar_rewrite_is_recorded_with_before_and_after(tmp_path):
    _clean(tmp_path)
    _write(tmp_path, "INCAR", "ALGO = All\nISMEAR = 1\nSIGMA = 0.2\n")
    audit = audit_vasp_inputs(tmp_path)
    assert audit["modified_files"] == ["INCAR"]
    assert audit["changes"]["INCAR"]["ALGO"] == {"original": "Normal", "current": "All"}
    assert audit["warnings"] == ["custodian modified INCAR"]


def test_a_kpoints_rewrite_is_recorded(tmp_path):
    _clean(tmp_path)
    _write(tmp_path, "KPOINTS", KPOINTS.replace("8 8 8", "16 16 16"))
    audit = audit_vasp_inputs(tmp_path)
    assert audit["modified_files"] == ["KPOINTS"]
    assert "kpoints" in audit["changes"]["KPOINTS"]


def test_a_poscar_atom_swap_is_recorded(tmp_path):
    """The one correction that would mis-index every force constant."""

    _clean(tmp_path)
    swapped = POSCAR.replace("Al Ta\n1 1\n", "Ta Al\n1 1\n")
    _write(tmp_path, "POSCAR", swapped)
    audit = audit_vasp_inputs(tmp_path)
    assert audit["modified_files"] == ["POSCAR"]
    assert audit["changes"]["POSCAR"]["species_order"] == {
        "original": ["Al", "Ta"], "current": ["Ta", "Al"]
    }


def test_a_poscar_coordinate_nudge_is_recorded(tmp_path):
    _clean(tmp_path)
    nudged = POSCAR.replace("0.5 0.5 0.5", "0.501 0.5 0.5")
    _write(tmp_path, "POSCAR", nudged)
    audit = audit_vasp_inputs(tmp_path)
    assert audit["modified_files"] == ["POSCAR"]
    assert audit["changes"]["POSCAR"]["frac_coords_max_abs_diff"] > 1e-6


def test_a_potcar_change_is_recorded_by_hash_only(tmp_path):
    _clean(tmp_path)
    _write(tmp_path, "POTCAR", "PAW_PBE Al\nPAW_PBE Ta_pv\n")
    audit = audit_vasp_inputs(tmp_path)
    assert audit["modified_files"] == ["POTCAR"]
    assert audit["hashes"]["POTCAR"]["original"] != audit["hashes"]["POTCAR"]["current"]


def test_missing_originals_are_reported_not_failed(tmp_path):
    """A trimmed job directory cannot prove the inputs were untouched.

    That is a provenance gap, not a reason to throw away a completed matrix --
    the same policy as a custodian correction.
    """

    _write(tmp_path, "INCAR", "ALGO = Normal\n")
    _write(tmp_path, "KPOINTS", KPOINTS)
    _write(tmp_path, "POSCAR", POSCAR)
    audit = audit_vasp_inputs(tmp_path)
    assert audit["original_input_available"] is False
    assert set(audit["missing_original_files"]) == {
        "INCAR", "KPOINTS", "POSCAR", "POTCAR"
    }
    assert audit["modified_files"] == []


def test_gzipped_inputs_are_audited_like_plain_ones(tmp_path):
    """atomate2 gzips the job directory at the end of every run."""

    _clean(tmp_path, gz=True)
    audit = audit_vasp_inputs(tmp_path)
    assert audit["original_files"]["INCAR"] == "INCAR.orig.gz"
    assert audit["current_files"]["INCAR"] == "INCAR.gz"
    assert audit["modified_files"] == []

    _write(tmp_path, "INCAR", "ALGO = All\nISMEAR = 1\nSIGMA = 0.2\n", gz=True)
    audit = audit_vasp_inputs(tmp_path)
    assert audit["modified_files"] == ["INCAR"]
    assert audit["changes"]["INCAR"]["ALGO"]["current"] == "All"


def test_a_periodic_image_is_not_a_modification(tmp_path):
    """0.0 and 1.0 are the same site; a re-wrap is not a custodian rewrite."""

    _clean(tmp_path)
    wrapped = POSCAR.replace("0.5 0.5 0.5", "1.5 1.5 1.5")
    _write(tmp_path, "POSCAR", wrapped)
    audit = audit_vasp_inputs(tmp_path)
    assert audit["modified_files"] == []


def test_a_truncated_gzip_is_recorded_not_raised(tmp_path):
    """An interrupted gzip_output_folder must not kill the collector.

    The audit is provenance: it may report that it could not read a file, but
    it may not take down a completed DFPT run three jobs downstream.
    """

    _clean(tmp_path, gz=True)
    path = tmp_path / "INCAR.orig.gz"
    data = path.read_bytes()
    path.write_bytes(data[: len(data) // 2])
    audit = audit_vasp_inputs(tmp_path)
    assert audit["modified_files"] == ["INCAR"]
    assert "error" in audit["changes"]["INCAR"]
