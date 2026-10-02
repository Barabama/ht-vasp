"""Canonical TB2J magnetic contract and derived observables.

The module deliberately keeps the producer contract independent from jobflow.
It can therefore canonicalize a parsed ``exchange.out`` result as well as a
JSON result reconstructed by a downstream reader.  The top-level legacy
fields are retained; the versioned contract lives under ``magnetic_contract``.
"""

from __future__ import annotations

import copy
import hashlib
import math
from pathlib import Path
from typing import Any


MAGNETIC_CONTRACT_SCHEMA_VERSION = "htvasp-tb2j-magnetic/v2"
LEGACY_SCHEMA_VERSION = "emtdb-magnetic-dataset/v1"
J_UNIT = "meV"
MOMENT_UNIT = "mu_B"
HAMILTONIAN = "H=-sum(J_ij S_i dot S_j)"
PAIR_COUNTING = "outgoing_per_site"
SPIN_NORMALIZATION = "unit"
REVERSE_POLICY = "allow_single_outgoing_validate_reverse_if_present"
K_B_MEV_PER_K = 0.08617333262

_STATUSES = {"computed", "nonmagnetic", "not_applicable", "unsupported", "needs_input"}
_TRANSITIONS = {"auto", "Curie", "Neel"}
_DENOMINATORS = {"structure_sites", "all_moment_sites"}
_TOL = 1.0e-10


def empty_magnetic_contract() -> dict[str, Any]:
    """Return a valid, input-free contract skeleton."""

    return {
        "schema_version": MAGNETIC_CONTRACT_SCHEMA_VERSION,
        "calphad": {
            "status": None,
            "production_ready": False,
            "phase": None,
            "type_character": None,
            "sublattice_ratios": [],
            "site_to_sublattice": {},
            "endmember_components": [],
            "order_num": None,
            "structure_role": None,
        },
        "exchange_convention": {
            "Hamiltonian": HAMILTONIAN,
            "J_unit": J_UNIT,
            "positive_J": "ferromagnetic",
            "pair_counting": PAIR_COUNTING,
            "spin_normalization": SPIN_NORMALIZATION,
            "reverse_policy": REVERSE_POLICY,
            # These aliases preserve the vocabulary of the v1 producer fields.
            "hamiltonian": HAMILTONIAN,
            "j_unit": J_UNIT,
            "j_positive_is_ferromagnetic": True,
            "reverse_pair_policy": REVERSE_POLICY,
        },
        "tc_model": {
            "method": "mfa_sublattice",
            "transition_type": "auto",
            "cutoff": "all_pairs",
            "k_B": K_B_MEV_PER_K,
        },
        "bmagn": {
            "source": "atom_moments.w_magmom",
            "unit": MOMENT_UNIT,
            "aggregation": "mean_abs",
            "denominator": "all_moment_sites",
        },
        "provenance": {
            "producer": "ht-vasp",
            "contract": MAGNETIC_CONTRACT_SCHEMA_VERSION,
            "physical_approval_requires_context": True,
        },
    }


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _as_dict(value: Any) -> dict[str, Any]:
    return dict(value) if isinstance(value, dict) else {}


def _context_root(magnetic_context: Any) -> tuple[dict[str, Any], bool]:
    if not isinstance(magnetic_context, dict) or not magnetic_context:
        return {}, False
    root = dict(magnetic_context)
    nested = root.get("magnetic_contract")
    if isinstance(nested, dict):
        root = nested
    return root, True


def _make_contract(magnetic_context: Any) -> tuple[dict[str, Any], bool]:
    root, supplied = _context_root(magnetic_context)
    contract = empty_magnetic_contract()

    calphad_input = _as_dict(root.get("calphad"))
    if not calphad_input:
        calphad_input = {
            key: root[key]
            for key in contract["calphad"]
            if key in root
        }
    if "sublattice_ratios" not in calphad_input and "ratios" in calphad_input:
        calphad_input["sublattice_ratios"] = calphad_input["ratios"]
    if "structure_role" not in calphad_input and root.get("structure_role") is not None:
        calphad_input["structure_role"] = root["structure_role"]
    calphad = contract["calphad"]
    for key in calphad:
        if key in calphad_input:
            calphad[key] = copy.deepcopy(calphad_input[key])

    exchange_input = _as_dict(root.get("exchange_convention"))
    exchange = contract["exchange_convention"]
    for key in exchange:
        if key in exchange_input:
            exchange[key] = copy.deepcopy(exchange_input[key])

    tc_input = _as_dict(root.get("tc_model"))
    tc = contract["tc_model"]
    for key in tc:
        if key in tc_input:
            tc[key] = copy.deepcopy(tc_input[key])
    # Accept the concise context spelling used by callers.
    for key in tc:
        if key in root and key not in tc_input:
            tc[key] = copy.deepcopy(root[key])

    bmagn_input = _as_dict(root.get("bmagn"))
    bmagn = contract["bmagn"]
    for key in bmagn:
        if key in bmagn_input:
            bmagn[key] = copy.deepcopy(bmagn_input[key])
    for key in bmagn:
        if key in root and key not in bmagn_input:
            bmagn[key] = copy.deepcopy(root[key])

    provenance = contract["provenance"]
    provenance["context_provided"] = supplied
    if "provenance" in root and isinstance(root["provenance"], dict):
        provenance.update(copy.deepcopy(root["provenance"]))
    return contract, supplied


def _site_element(site: dict[str, Any]) -> str | None:
    species = site.get("species")
    if isinstance(species, list) and species and isinstance(species[0], dict):
        value = species[0].get("element") or species[0].get("label")
        if isinstance(value, str) and value.strip():
            return value.strip().split(" ", 1)[0]
    label = site.get("label")
    if isinstance(label, str) and label.strip():
        letters = ""
        for char in label.strip():
            if not char.isalpha():
                break
            letters += char
        return letters or None
    return None


def _moment_label(moment: dict[str, Any], index: int) -> str:
    for key in ("atom", "label", "site"):
        value = moment.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return str(index)


def _moment_value(moment: dict[str, Any]) -> float | None:
    for key in ("w_magmom", "magmom", "moment"):
        value = moment.get(key)
        if value is not None:
            try:
                return float(value)
            except (TypeError, ValueError):
                return None
    return None


def _moments(data: dict[str, Any]) -> tuple[list[tuple[str, float]], list[str]]:
    values: list[tuple[str, float]] = []
    missing: list[str] = []
    raw = data.get("atom_moments")
    if not isinstance(raw, list):
        return values, ["atom_moments"]
    for index, item in enumerate(raw):
        if not isinstance(item, dict):
            missing.append(f"atom_moments[{index}]")
            continue
        value = _moment_value(item)
        if value is None:
            missing.append(f"atom_moments[{index}].w_magmom")
            continue
        values.append((_moment_label(item, index), value))
    return values, missing


def _structure_site_labels(data: dict[str, Any], moment_values: list[tuple[str, float]]) -> list[str]:
    structure = data.get("structure")
    sites = structure.get("sites") if isinstance(structure, dict) else None
    labels: list[str] = []
    if isinstance(sites, list):
        element_counts: dict[str, int] = {}
        for index, site in enumerate(sites):
            if not isinstance(site, dict):
                labels.append(str(index))
                continue
            explicit = site.get("label")
            if isinstance(explicit, str) and explicit.strip():
                labels.append(explicit.strip())
                continue
            element = _site_element(site)
            if element:
                element_counts[element] = element_counts.get(element, 0) + 1
                labels.append(f"{element}{element_counts[element]}")
            else:
                labels.append(str(index))
    if not labels:
        labels = [label for label, _ in moment_values]
    return labels


def _ratio_names(calphad: dict[str, Any], site_map: Any) -> list[str]:
    ratios = calphad.get("sublattice_ratios")
    if isinstance(ratios, dict):
        names = [str(key) for key in ratios]
    elif isinstance(ratios, list):
        explicit = calphad.get("sublattice_names")
        if isinstance(explicit, list) and len(explicit) == len(ratios):
            names = [str(value) for value in explicit]
        elif isinstance(site_map, dict):
            names = []
            for value in site_map.values():
                if isinstance(value, (int, float)) and float(value).is_integer():
                    candidate = str(int(value))
                else:
                    candidate = str(value)
                if candidate not in names:
                    names.append(candidate)
            names = names[: len(ratios)]
        else:
            names = [str(index) for index in range(len(ratios))]
    else:
        names = []
    return names


def _site_map(calphad: dict[str, Any], data: dict[str, Any], moment_values: list[tuple[str, float]]) -> tuple[dict[str, str], list[str], list[str]]:
    raw_map = calphad.get("site_to_sublattice")
    structure_labels = _structure_site_labels(data, moment_values)
    moment_labels = [label for label, _ in moment_values]
    labels = list(dict.fromkeys(structure_labels + moment_labels))
    if isinstance(raw_map, list):
        raw_map = {str(index): value for index, value in enumerate(raw_map)}
    if not isinstance(raw_map, dict):
        return {}, [], labels

    names = _ratio_names(calphad, raw_map)
    if not names:
        for value in raw_map.values():
            name = str(value)
            if name not in names:
                names.append(name)

    result: dict[str, str] = {}

    def normalize(value: Any) -> str:
        if isinstance(value, (int, float)) and float(value).is_integer():
            index = int(value)
            return names[index] if 0 <= index < len(names) else str(index)
        text = str(value)
        if text.isdigit() and int(text) < len(names):
            return names[int(text)]
        return text

    by_index = {str(index): label for index, label in enumerate(labels)}
    for key, value in raw_map.items():
        key_text = str(key)
        label = by_index.get(key_text, key_text)
        result[label] = normalize(value)

    # A map keyed by elements is useful for one-site-per-element endmembers.
    for label in labels:
        if label in result:
            continue
        element = "".join(char for char in label if char.isalpha())
        if element in result:
            result[label] = result[element]
    names = list(dict.fromkeys(names + list(result.values())))
    return result, names, labels


def _cutoff_value(value: Any) -> float | None:
    if value is None or value == "all_pairs" or value == "all":
        return None
    try:
        cutoff = float(value)
    except (TypeError, ValueError):
        return None
    return cutoff if math.isfinite(cutoff) and cutoff >= 0 else None


def _pair_j(pair: dict[str, Any]) -> float | None:
    for key in ("J_iso_meV", "J", "j"):
        if key in pair:
            try:
                value = float(pair[key])
            except (TypeError, ValueError):
                return None
            return value if math.isfinite(value) else None
    return None


def _pair_distance(pair: dict[str, Any]) -> float | None:
    try:
        return float(pair["distance"])
    except (KeyError, TypeError, ValueError):
        return None


def _anisotropy_warnings(data: dict[str, Any]) -> list[str]:
    pairs = data.get("j_pairs") if isinstance(data.get("j_pairs"), list) else []
    tensors = data.get("j_tensors") if isinstance(data.get("j_tensors"), list) else []
    warnings: list[str] = []
    for index, tensor in enumerate(tensors[: len(pairs)]):
        if not isinstance(tensor, list) or len(tensor) != 3 or any(
            not isinstance(row, list) or len(row) != 3 for row in tensor
        ):
            continue
        scalar = _pair_j(pairs[index]) if isinstance(pairs[index], dict) else None
        if scalar is None:
            continue
        tolerance = max(1.0e-3, abs(scalar) * 1.0e-3)
        anisotropic = any(
            abs(float(tensor[row][column]) - (scalar if row == column else 0.0)) > tolerance
            for row in range(3)
            for column in range(3)
        )
        if anisotropic:
            warnings.append(
                f"anisotropic exchange tensor at pair index {index}; MFA uses J_iso_meV"
            )
    return warnings


def _jacobi_eigenpairs(matrix: list[list[float]]) -> tuple[list[float], list[float]] | None:
    """Return eigenvalues and the eigenvector for the largest one.

    The exchange matrix is expected to be symmetric after the per-site
    averaging.  Jacobi keeps this module dependency-free and is sufficient for
    the small sublattice matrices used by the producer contract.
    """

    size = len(matrix)
    if size == 0:
        return [], []
    if size == 1:
        return [float(matrix[0][0])], [1.0]
    if any(abs(matrix[i][j] - matrix[j][i]) > 1.0e-7 for i in range(size) for j in range(size)):
        return None

    values = [[float(value) for value in row] for row in matrix]
    vectors = [[1.0 if row == column else 0.0 for column in range(size)] for row in range(size)]
    max_iterations = max(32, 32 * size * size)
    for _ in range(max_iterations):
        pivot = max(
            ((abs(values[i][j]), i, j) for i in range(size) for j in range(i + 1, size)),
            default=(0.0, 0, 0),
        )
        magnitude, p, q = pivot
        if magnitude <= 1.0e-12:
            break
        angle = 0.5 * math.atan2(2.0 * values[p][q], values[q][q] - values[p][p])
        cosine, sine = math.cos(angle), math.sin(angle)
        for index in range(size):
            if index not in (p, q):
                vip, viq = values[index][p], values[index][q]
                values[index][p] = values[p][index] = cosine * vip - sine * viq
                values[index][q] = values[q][index] = sine * vip + cosine * viq
        app, aqq, apq = values[p][p], values[q][q], values[p][q]
        values[p][p] = cosine * cosine * app - 2.0 * sine * cosine * apq + sine * sine * aqq
        values[q][q] = sine * sine * app + 2.0 * sine * cosine * apq + cosine * cosine * aqq
        values[p][q] = values[q][p] = 0.0
        for index in range(size):
            vip, viq = vectors[index][p], vectors[index][q]
            vectors[index][p] = cosine * vip - sine * viq
            vectors[index][q] = sine * vip + cosine * viq
    eigenvalues = [values[index][index] for index in range(size)]
    largest = max(range(size), key=lambda index: eigenvalues[index])
    vector = [vectors[index][largest] for index in range(size)]
    return eigenvalues, vector


def _two_by_two_eigenpairs(matrix: list[list[float]]) -> tuple[list[float], list[float]] | None:
    a, b = matrix[0]
    c, d = matrix[1]
    discriminant = (a - d) ** 2 + 4.0 * b * c
    if discriminant < -1.0e-9:
        return None
    root = math.sqrt(max(0.0, discriminant))
    eigenvalues = [(a + d - root) / 2.0, (a + d + root) / 2.0]
    largest = 1
    eigenvalue = eigenvalues[largest]
    if abs(b) > _TOL:
        vector = [b, eigenvalue - a]
    elif abs(c) > _TOL:
        vector = [eigenvalue - d, c]
    else:
        vector = [1.0, 0.0]
    return eigenvalues, vector


def _largest_eigenpair(matrix: list[list[float]]) -> tuple[list[float], float, list[float]] | None:
    if len(matrix) == 2:
        result = _two_by_two_eigenpairs(matrix)
    else:
        result = _jacobi_eigenpairs(matrix)
    if result is None:
        return None
    eigenvalues, vector = result
    index = max(range(len(eigenvalues)), key=lambda item: eigenvalues[item])
    if len(matrix) == 2 and index == 0:
        # Recompute the eigenvector for the smaller root only when necessary.
        a, b = matrix[0]
        c, d = matrix[1]
        eigenvalue = eigenvalues[index]
        if abs(b) > _TOL:
            vector = [b, eigenvalue - a]
        elif abs(c) > _TOL:
            vector = [eigenvalue - d, c]
        else:
            vector = [1.0, 0.0]
    return eigenvalues, eigenvalues[index], vector


def _transition_from_vector(vector: list[float]) -> str | None:
    signs = [value > 1.0e-8 for value in vector if abs(value) > 1.0e-8]
    negative = [value < -1.0e-8 for value in vector if abs(value) > 1.0e-8]
    if len(signs) != len(vector) or not signs:
        return None
    if any(signs) and any(negative):
        return "Neel"
    return "Curie"


def _bmagn(
    data: dict[str, Any],
    contract: dict[str, Any],
    moment_values: list[tuple[str, float]],
) -> tuple[float | None, str | None]:
    bmagn = contract["bmagn"]
    denominator = bmagn.get("denominator")
    aggregation = bmagn.get("aggregation")
    if aggregation not in {"mean_abs", "mean_abs_w_magmom", "mean_abs_per_atom"}:
        return None, f"unsupported BMAGN aggregation: {aggregation!r}"
    if denominator not in _DENOMINATORS:
        return None, f"unsupported BMAGN denominator: {denominator!r}"
    if not moment_values:
        return None, "no local magnetic moments"
    denominator_count = len(moment_values)
    if denominator == "structure_sites":
        structure = data.get("structure")
        sites = structure.get("sites") if isinstance(structure, dict) else None
        denominator_count = len(sites) if isinstance(sites, list) else int(data.get("num_atoms") or 0)
    if denominator_count <= 0:
        return None, "BMAGN denominator has no sites"
    return sum(abs(value) for _, value in moment_values) / denominator_count, None


def derive_magnetic_observables(
    data: dict[str, Any],
    contract: dict[str, Any],
    *,
    context_supplied: bool,
    warnings: list[str],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Derive MFA sublattice observables and return ``(observables, audit)``."""

    moment_values, moment_missing = _moments(data)
    bmagn_value, bmagn_error = _bmagn(data, contract, moment_values)
    calphad = contract["calphad"]
    tc_model = contract["tc_model"]
    site_map, sublattices, all_site_labels = _site_map(calphad, data, moment_values)
    pairs = data.get("j_pairs") if isinstance(data.get("j_pairs"), list) else []
    cutoff = _cutoff_value(tc_model.get("cutoff"))
    selected_pairs: list[dict[str, Any]] = []
    invalid_pairs = 0
    excluded_pairs = 0
    for pair in pairs:
        if not isinstance(pair, dict) or _pair_j(pair) is None:
            invalid_pairs += 1
            continue
        distance = _pair_distance(pair)
        if cutoff is not None and (distance is None or distance > cutoff):
            excluded_pairs += 1
            continue
        selected_pairs.append(pair)

    source_labels = list(dict.fromkeys(
        [str(pair.get("i")) for pair in selected_pairs if pair.get("i") is not None]
        + all_site_labels
    ))
    missing_sites = sorted({
        str(pair.get(key))
        for pair in selected_pairs
        for key in ("i", "j")
        if pair.get(key) is not None and str(pair.get(key)) not in site_map
    })

    site_sums: dict[tuple[str, str], float] = {}
    for pair in selected_pairs:
        source, target = str(pair.get("i")), str(pair.get("j"))
        source_sublattice = site_map.get(source)
        target_sublattice = site_map.get(target)
        if source_sublattice is None or target_sublattice is None:
            continue
        key = (source, target_sublattice)
        site_sums[key] = site_sums.get(key, 0.0) + float(_pair_j(pair))

    sites_by_sublattice: dict[str, list[str]] = {name: [] for name in sublattices}
    for label in list(dict.fromkeys(source_labels + list(site_map))):
        sublattice = site_map.get(label)
        if sublattice in sites_by_sublattice and label not in sites_by_sublattice[sublattice]:
            sites_by_sublattice[sublattice].append(label)

    matrix: list[list[float]] = []
    for source_sublattice in sublattices:
        source_sites = sites_by_sublattice.get(source_sublattice, [])
        row: list[float] = []
        for target_sublattice in sublattices:
            if not source_sites:
                row.append(0.0)
                continue
            row.append(
                sum(site_sums.get((label, target_sublattice), 0.0) for label in source_sites)
                / len(source_sites)
            )
        matrix.append(row)

    trace: dict[str, Any] = {
        "formula": "TC_or_TN=(2/3)*lambda_max/k_B",
        "k_B": float(tc_model.get("k_B", K_B_MEV_PER_K)),
        "k_B_unit": "meV/K",
        "J_unit": J_UNIT,
        "pair_count_total": len(pairs),
        "pair_count_used": len(selected_pairs),
        "pair_count_excluded_by_cutoff": excluded_pairs,
        "invalid_pair_count": invalid_pairs,
        "cutoff": tc_model.get("cutoff"),
        "pair_counting": PAIR_COUNTING,
        "matrix_construction": (
            "sum outgoing J per source site to each target sublattice, then average source sites"
        ),
        "sublattices": list(sublattices),
        "sublattice_site_counts": {
            name: len(sites_by_sublattice.get(name, [])) for name in sublattices
        },
        "source_site_count": len(source_labels),
        "tensor_policy": "J_iso_meV is used; anisotropic tensors are retained and warned",
    }
    if missing_sites:
        trace["unmapped_sites"] = missing_sites
    if moment_missing:
        trace["moment_input_issues"] = moment_missing
    if bmagn_error:
        trace["bmagn_issue"] = bmagn_error

    observable: dict[str, Any] = {
        "status": "needs_input",
        "TC": None,
        "TN": None,
        "BMAGN": bmagn_value,
        "exchange_matrix": matrix if sublattices else [],
        "exchange_matrix_unit": J_UNIT,
        "eigenvalues": [],
        "transition_type": None,
        "conversion_trace": trace,
    }
    missing_context: list[str] = []
    for key in ("phase", "type_character", "sublattice_ratios", "site_to_sublattice", "endmember_components", "order_num", "structure_role"):
        value = calphad.get(key)
        if value is None or value == "" or value == [] or value == {}:
            missing_context.append(f"calphad.{key}")
    if tc_model.get("method") != "mfa_sublattice":
        missing_context.append("tc_model.method=mfa_sublattice")
    if tc_model.get("transition_type") not in _TRANSITIONS:
        missing_context.append("tc_model.transition_type")
    if not isinstance(tc_model.get("k_B"), (int, float)) or float(tc_model["k_B"]) <= 0:
        missing_context.append("tc_model.k_B")
    if not sublattices:
        missing_context.append("calphad.sublattice_ratios")
    if missing_sites:
        missing_context.append("calphad.site_to_sublattice for all J sites")
    if bmagn_error and "no local magnetic moments" not in bmagn_error:
        missing_context.append("bmagn")
    if not moment_values:
        missing_context.append("atom_moments.w_magmom")

    all_moments_zero = bool(moment_values) and all(abs(value) <= _TOL for _, value in moment_values)
    if context_supplied and calphad.get("status") == "not_applicable":
        observable["status"] = "not_applicable"
        trace["not_applicable_reason"] = (
            "producer context marks this structure as a reference state; "
            "magnetic CALPHAD parameters are supplied by the unary database"
        )
    elif not context_supplied:
        observable["status"] = "needs_input"
        trace["nonproduction_reason"] = "magnetic_context is absent"
    elif missing_context:
        observable["status"] = "needs_input"
        trace["missing_context"] = missing_context
    elif all_moments_zero:
        observable["status"] = "nonmagnetic"
        trace["nonmagnetic_reason"] = "all local moments are zero"
    else:
        eigen = _largest_eigenpair(matrix)
        if eigen is None:
            observable["status"] = "unsupported"
            trace["unsupported_reason"] = "exchange matrix has no real eigen-decomposition"
        else:
            eigenvalues, lambda_max, vector = eigen
            observable["eigenvalues"] = [float(value) for value in eigenvalues]
            trace["lambda_max_meV"] = float(lambda_max)
            transition = _transition_from_vector(vector) if lambda_max > _TOL else None
            if lambda_max <= _TOL:
                observable["status"] = "nonmagnetic"
                trace["nonmagnetic_reason"] = "maximum real eigenvalue is non-positive"
            elif transition is None:
                observable["status"] = "unsupported"
                trace["unsupported_reason"] = "largest-eigenvalue vector has no uniform sign pattern"
            else:
                temperature = (2.0 / 3.0) * lambda_max / float(tc_model["k_B"])
                observable["status"] = "computed"
                observable["transition_type"] = transition
                if transition == "Curie":
                    observable["TC"] = temperature
                else:
                    observable["TN"] = temperature
    if observable["status"] == "needs_input" and not context_supplied:
        observable["eigenvalues"] = []
    if observable["status"] in {"computed", "nonmagnetic"} and not observable["eigenvalues"] and matrix:
        eigen = _largest_eigenpair(matrix)
        if eigen is not None:
            observable["eigenvalues"] = [float(value) for value in eigen[0]]

    audit = {
        "contract_status": (
            "raw_reader_only" if not context_supplied else observable["status"]
        ),
        "production_ready": bool(
            context_supplied and observable["status"] in {"computed", "nonmagnetic", "not_applicable"}
        ),
        "review_required": bool(
            not context_supplied
            or observable["status"] not in {"computed", "nonmagnetic", "not_applicable"}
            or warnings
        ),
        "verified_checks": [
            "magnetic contract schema version recorded",
            "H=-sum(J_ij S_i dot S_j), J in meV, unit spins recorded",
            "all source J records retained; cutoff only affects derived matrix",
            "raw reader contract is not physical production approval",
        ],
    }
    if bmagn_value is not None:
        audit["verified_checks"].append(
            f"BMAGN=mean(abs(w_magmom)) with denominator={contract['bmagn']['denominator']}"
        )
    if matrix and not missing_sites:
        audit["verified_checks"].append("sublattice exchange matrix averaged per source site")
    if warnings:
        audit["verified_checks"].append("anisotropic tensor warning recorded")
    if trace.get("lambda_max_meV") is not None:
        audit["verified_checks"].append("MFA temperature conversion and lambda_max recorded")
    return observable, audit


def _source_metadata(result: dict[str, Any], output_path: str | Path | None) -> tuple[Path | None, dict[str, Any]]:
    path_value = result.get("exchange_out_path")
    source_path = Path(str(path_value)).resolve() if path_value else None
    if output_path is not None:
        candidate = Path(output_path).resolve().parent / "tb2jflow" / "3-tb2j_solve" / "exchange.out"
        if candidate.is_file():
            source_path = candidate
    source = _as_dict(result.get("source"))
    if source_path is not None:
        result["exchange_out_path"] = source_path.as_posix()
        source.update({"path": source_path.as_posix(), "format": "exchange_out"})
        if source_path.is_file():
            source["sha256"] = _sha256_file(source_path)
    source.setdefault("format", "exchange_out")
    source.setdefault("producer", result.get("producer", "TB2J"))
    source.setdefault("producer_version", result.get("producer_version"))
    return source_path, source


def canonicalize_tb2j_result(
    data: dict[str, Any],
    output_path: str | Path | None = None,
    magnetic_context: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a v2 magnetic-contract JSON-compatible TB2J result.

    ``magnetic_context=None`` is intentionally valid.  It produces a raw
    canonical artifact with ``needs_input`` and ``production_ready=False``;
    parser provenance alone never grants physical production approval.
    """

    if not isinstance(data, dict):
        raise ValueError("TB2J result must be a dictionary")
    result = copy.deepcopy(data)
    result.setdefault("schema_version", LEGACY_SCHEMA_VERSION)
    pairs = result.get("j_pairs")
    if not isinstance(pairs, list):
        pairs = []
        result["j_pairs"] = pairs
    tensors = result.get("j_tensors")
    if not isinstance(tensors, list):
        tensors = []
        result["j_tensors"] = tensors
    if tensors and len(tensors) != len(pairs):
        raise ValueError("TB2J result requires j_tensors parallel to j_pairs")
    result["num_j_pairs"] = len(pairs)
    moments = result.get("atom_moments")
    if not isinstance(moments, list):
        result["atom_moments"] = []
        moments = result["atom_moments"]
    if not result.get("num_atoms") and moments:
        result["num_atoms"] = len(moments)

    structure = result.get("structure")
    if isinstance(structure, dict) and isinstance(structure.get("sites"), list):
        sites = structure["sites"]
        site_records: list[tuple[str | None, str | None]] = []
        by_label: dict[str, str] = {}
        for index, site in enumerate(sites):
            if not isinstance(site, dict):
                site_records.append((None, None))
                continue
            site_element = _site_element(site)
            site_label = site.get("label") or site.get("id")
            normalized_label = str(site_label).strip() if site_label else None
            site_records.append((normalized_label, site_element))
            if normalized_label and site_element:
                by_label[normalized_label] = site_element

        # TB2J's Atoms section and the POSCAR/Structure site order are both
        # valid orderings, but they are not required to be identical.  Match
        # by the stable atom label first; only unlabeled records use positional
        # fallback.
        for index, (site_label, positional_element) in enumerate(site_records):
            moment = moments[index] if index < len(moments) else None
            if not isinstance(moment, dict):
                continue
            moment_label = _moment_label(moment, index)
            moment_element = "".join(char for char in moment_label if char.isalpha())
            site_element = by_label.get(moment_label)
            if site_element is None and site_label is None:
                site_element = positional_element
            if site_element and moment_element and site_element.lower() != moment_element.lower():
                raise ValueError(
                    f"TB2J structure site {index} element {site_element!r} does not match "
                    f"atom_moments element {moment_element!r}"
                )

        # The native TB2J ``Atoms`` section can use a different order from
        # POSCAR.  Canonical JSON keeps the explicit atom-moment order so the
        # structure/moment arrays are positionally aligned for downstream
        # readers, while labels preserve the original site identity.
        site_by_label = {
            str(site.get("label") or site.get("id")): site
            for site in sites
            if isinstance(site, dict) and (site.get("label") or site.get("id"))
        }
        moment_labels = [_moment_label(moment, index) for index, moment in enumerate(moments)]
        if len(site_by_label) == len(sites) and set(moment_labels) == set(site_by_label):
            structure["sites"] = [site_by_label[label] for label in moment_labels]

    source_path, source = _source_metadata(result, output_path)
    result["source"] = source
    provenance = _as_dict(result.get("provenance"))
    provenance.update({
        "source_format": source.get("format"),
        "source_path": source.get("path"),
        "source_sha256": source.get("sha256"),
        "producer": source.get("producer"),
        "producer_version": source.get("producer_version"),
        "j_unit": J_UNIT,
        "moment_unit": MOMENT_UNIT,
        "hamiltonian": HAMILTONIAN,
        "j_positive_is_ferromagnetic": True,
        "pair_counting_rule": PAIR_COUNTING,
        "spin_normalization": SPIN_NORMALIZATION,
        "reverse_pair_policy": REVERSE_POLICY,
        "magnetic_contract": MAGNETIC_CONTRACT_SCHEMA_VERSION,
    })

    contract, context_supplied = _make_contract(magnetic_context)
    warnings = list(result.get("warnings") or [])
    for warning in _anisotropy_warnings(result):
        if warning not in warnings:
            warnings.append(warning)
    root, _ = _context_root(magnetic_context)
    fail_on_anisotropy = bool(root.get("fail_on_anisotropy"))
    if isinstance(root.get("exchange_convention"), dict):
        fail_on_anisotropy = bool(root["exchange_convention"].get("fail_on_anisotropy", fail_on_anisotropy))
    if warnings and fail_on_anisotropy and any("anisotropic exchange tensor" in warning for warning in warnings):
        raise ValueError("anisotropic exchange tensor rejected by fail_on_anisotropy")

    observables, audit = derive_magnetic_observables(
        result,
        contract,
        context_supplied=context_supplied,
        warnings=warnings,
    )
    audit["source_path"] = source.get("path")
    audit["source_sha256"] = source.get("sha256")
    audit["parser_audit"] = _as_dict(result.get("audit"))
    contract["provenance"].update({
        "source_format": source.get("format"),
        "source_path": source.get("path"),
        "source_sha256": source.get("sha256"),
        "context_provided": context_supplied,
        "derived_observables": ["TC", "TN", "BMAGN", "exchange_matrix", "eigenvalues"],
    })
    result.update({
        "warnings": warnings,
        "magnetic_contract": contract,
        "magnetic_observables": observables,
        "audit": audit,
        "provenance": provenance,
    })
    for moment in moments:
        if isinstance(moment, dict):
            moment.setdefault("unit", MOMENT_UNIT)
            moment.setdefault("moment_unit", MOMENT_UNIT)
    return result


__all__ = [
    "HAMILTONIAN",
    "J_UNIT",
    "K_B_MEV_PER_K",
    "MAGNETIC_CONTRACT_SCHEMA_VERSION",
    "MOMENT_UNIT",
    "PAIR_COUNTING",
    "REVERSE_POLICY",
    "canonicalize_tb2j_result",
    "derive_magnetic_observables",
    "empty_magnetic_contract",
]
