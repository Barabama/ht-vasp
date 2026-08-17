"""TB2J Input Set Generator - VASP inputs for SCF + Wannier90 magnetic exchange."""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from atomate2.vasp.sets.base import VaspInputGenerator
from pymatgen.io.vasp import Kpoints
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from pymatgen.core import Structure
    from pymatgen.io.vasp import VaspInput

log = logging.getLogger(__name__)

# Number of Wannier orbitals assumed per atom (s + p + d = 1 + 3 + 5).
# Validated for 3d magnetic systems: NUM_WANN = n_atoms * 9.
ORBITALS_PER_ATOM = 9

# Wannier90-interface INCAR settings (from verified submit_tb2j.sh INCAR_W90).
# These are merged into the INCAR on top of the inherited user_incar_settings,
# so a user may still override any of them. NUM_WANN / NBANDS are set
# dynamically in get_input_set().
W90_INCAR_DEFAULTS: dict[str, Any] = {
    "LWANNIER90": True,  # run the Wannier90 interface inside VASP
    "LWANNIER90_RUN": True,  # execute wannier90.x from VASP
    "LWRITE_MMN_AMN": True,  # write overlap/projection files for TB2J
    "KPAR": 1,  # Wannier90 interface requires KPAR = 1
    "LREAL": False,  # reciprocal-space projections
    "NSW": 0,  # static SCF only
    "IBRION": -1,  # static SCF only
    "NELM": 200,
    "LORBIT": 11,
}


@dataclass
class Tb2jInputSetGenerator(VaspInputGenerator):
    """Generator for SCF + Wannier90 input sets feeding a TB2J exchange run.

    Inherits from :class:`VaspInputGenerator` for atomate2 compatibility. The
    ``get_input_set`` override injects the Wannier90-interface INCAR keys
    (LWANNIER90, KPAR=1, NUM_WANN, NBANDS, ...) and, as a side effect, writes a
    ``wannier90.win`` file into the job directory so VASP's Wannier90 interface
    picks up the explicit orbital projections (VASP's auto-generated win lacks
    ``begin projections`` and fails with "No projection mode").

    Args:
        elements: Magnetic elements used for the TB2J ``--elements`` flag and,
            when no structure is available, as the fallback projection list.
        num_wann: Number of Wannier orbitals. If None, auto-derived as
            ``n_atoms * ORBITALS_PER_ATOM`` (9 orbitals per spd atom).
        num_bands_extra: Extra bands above num_wann before the NPAR rounding
            (``nbands_min = num_wann + this``, see ``get_num_bands``).
        ntasks: Number of MPI ranks the SCF runs on. With KPAR=1 and NCORE=1
            (forced by the Wannier90 interface) VASP uses ``NPAR = NTASKS`` and
            checks that NBANDS is divisible by NPAR, silently bumping it to the
            next multiple otherwise. NBANDS is therefore pre-rounded to a
            multiple of ``ntasks`` so the INCAR NBANDS always equals the
            wannier90.win ``num_bands`` (no ``The number of bands has been
            changed`` drift / win mismatch).
        kpoints_grid: Explicit k-mesh grid for the Wannier90 SCF. Wannier90
            requires an explicit, isotropic k-mesh (the atomate2 auto-kspacing /
            reciprocal-density mesh can be non-isotropic and fails with
            ``kmesh_get_bvector: Not enough bvectors found``). Defaults to the
            verified 8x8x8 Gamma-centered grid from ``tb2j_test/SER-Co/W90``.
    """

    elements: list[str] = field(default_factory=list)
    num_wann: int | None = None
    num_bands_extra: int = 20
    ntasks: int = 32
    # NOTE: not named `kpoints` - the base VaspInputSet exposes a read-only
    # `kpoints` property (generates the KPOINTS file via kpoints_updates). A
    # dataclass field with that name would shadow it with a plain tuple.
    kpoints_grid: tuple = (8, 8, 8)

    def __post_init__(self) -> None:
        super().__post_init__()

    @property
    def kpoints_updates(self) -> Kpoints:
        """Override the atomate2 auto k-mesh with an explicit Gamma grid.

        Wannier90 needs an explicit, isotropic k-mesh. The atomate2 default
        (``reciprocal_density`` based, ``auto_kspacing``) produces a symmetry-
        reduced mesh such as 10x10x6 for HCP SER-Co, which breaks the W90
        interface with ``kmesh_get_bvector: Not enough bvectors found``. Returning
        a ``Kpoints`` object from this property makes pymatgen use it directly
        (``kpoints_updates`` accepts a dict or a ``Kpoints`` object) and ensures
        KSPACING is unset in the INCAR (KPOINTS file takes precedence).

        Note: ignored if ``user_kpoints_settings`` is set explicitly.
        """
        return Kpoints.gamma_automatic(kpts=tuple(self.kpoints_grid))

    def get_num_wann(self, structure: Structure | None = None) -> int:
        """Resolve the number of Wannier orbitals.

        Priority: explicit ``num_wann`` > auto ``n_atoms * ORBITALS_PER_ATOM``.

        Args:
            structure: Structure used for auto-derivation.

        Returns:
            Number of Wannier orbitals.
        """
        if self.num_wann is not None:
            return self.num_wann
        if structure is None:
            structure = self.structure
        if structure is None:
            raise ValueError(
                "Cannot derive NUM_WANN: no structure and num_wann not set. "
                "Provide a structure or set num_wann explicitly."
            )
        return len(structure) * ORBITALS_PER_ATOM

    def get_num_bands(self, num_wann: int | None = None) -> int:
        """Resolve NBANDS rounded up to a multiple of NPAR (= NTASKS).

        VASP requires NBANDS to be divisible by NPAR. With the Wannier90
        interface (KPAR=1, NCORE=1) ``NPAR = NTASKS``, and VASP silently bumps
        an indivisible NBANDS up to the next multiple, printing ``number of
        bands changed``. Because wannier90.win ``num_bands`` is a fixed value,
        that drift makes the INCAR NBANDS and the win disagree. Pre-rounding
        here (``NBANDS = ceil((num_wann + num_bands_extra) / ntasks) * ntasks``)
        keeps VASP from changing bands, so the INCAR and the win always match.

        Args:
            num_wann: Number of Wannier orbitals. If None, resolved via
                ``get_num_wann`` (requires an explicit ``num_wann`` or a
                structure).

        Returns:
            NBANDS as the smallest multiple of ``ntasks`` >= num_wann + extra.
        """
        num_wann = self.get_num_wann() if num_wann is None else num_wann
        npar = self.ntasks  # KPAR=1, NCORE=1 -> NPAR = NTASKS
        nbands_min = num_wann + self.num_bands_extra
        return int(math.ceil(nbands_min / npar) * npar)

    def _projection_elements(self, structure: Structure | None = None) -> list[str]:
        """Get element symbols for the win projections in POSCAR site order.

        Args:
            structure: Structure to read element order from.

        Returns:
            List of element symbols (one entry per distinct element, site order).
        """
        if structure is not None:
            seen: list[str] = []
            for site in structure:
                sym = site.species_string
                if sym not in seen:
                    seen.append(sym)
            return seen
        if self.elements:
            return list(self.elements)
        raise ValueError(
            "Cannot build wannier90.win projections: provide a structure or "
            "set the `elements` list."
        )

    def _symmetrize_structure(self, structure: Structure) -> Structure:
        """Snap a near-symmetric lattice to the detected spacegroup symmetry.

        The preceding R3 relaxation can leave a tiny numerical asymmetry in
        the relaxed lattice (e.g. HCP SER-Co at a1.y = -3e-6, angle 120.0001
        deg instead of 120 deg). This is most pronounced when ISYM=0 is left
        in the global INCAR; the Tb2jWorker now forces ISYM=1 for the R3
        relax, but the guard is kept as a safety net for any residual
        numerical noise (and for callers that do not force ISYM=1).
        wannier90's k-mesh setup then splits otherwise-degenerate neighbour
        shells and fails with ``kmesh_get_bvector: Not enough bvectors
        found`` even for the explicit Gamma 8x8x8 grid. Refining to the
        detected spacegroup restores the shell degeneracy (verified: perfect
        hexagonal lattice passes wannier90 kmesh setup for the same grid).

        The refinement is applied only when it preserves the cell (same number
        of sites and same composition); otherwise the structure is used as-is
        so genuinely low-symmetry cells are not altered.

        The symmetry tolerance is escalated (1e-3 -> 2e-2) only when the tight
        tolerance changes the atom count: a tight tolerance can detect a
        spuriously low-symmetry group for a slightly distorted high-symmetry
        cell. E.g. an R3-relaxed BCC Fe cell whose angle gamma drifted from the
        ideal 109.47 deg to 109.89 deg is detected as Fmmm (69) at symprec=1e-3,
        whose conventional cell holds 4 atoms (so the atom-count guard rejects
        it and wannier90 still fails with ``kmesh_get_bvector``). At symprec=2e-2
        the same cell is detected as Im-3m (229), and ``get_refined_structure``
        returns the clean 2-atom cubic conventional BCC cell with all 8 NN
        shells exactly degenerate (verified). The first tolerance that preserves
        the atom count is accepted.

        Args:
            structure: Structure to symmetrize (typically the R3-relaxed one).

        Returns:
            Symmetrized structure, or the input unchanged if no refinement
            preserves the cell.
        """
        for symprec in (1e-3, 5e-3, 1e-2, 2e-2):
            try:
                refined = SpacegroupAnalyzer(
                    structure, symprec=symprec
                ).get_refined_structure()
            except Exception:
                log.debug(
                    f"Spacegroup refinement failed at symprec={symprec}",
                    exc_info=True,
                )
                continue
            if (
                len(refined) != len(structure)
                or refined.composition != structure.composition
            ):
                log.debug(
                    f"Symmetrization at symprec={symprec} would change the cell "
                    f"({len(structure)} -> {len(refined)} sites); trying looser "
                    f"tolerance"
                )
                continue
            if refined.lattice != structure.lattice:
                log.debug(
                    f"Symmetrized lattice (symprec={symprec}): "
                    f"{structure.lattice.abc} / {structure.lattice.angles} -> "
                    f"{refined.lattice.abc} / {refined.lattice.angles}"
                )
            return refined
        log.debug("No symmetry refinement preserved the cell; using structure as-is")
        return structure

    def get_wannier90_win_string(
        self,
        structure: Structure | None = None,
        num_wann: int | None = None,
    ) -> str:
        """Generate the ``wannier90.win`` content.

        Template validated on SER-Co / BCC-Co-Ni / FCC-Fe-Mn (submit_tb2j.sh).

        Args:
            structure: Structure used for projection order and num_wann derivation.
            num_wann: Explicit num_wann override (else resolved via get_num_wann).

        Returns:
            wannier90.win file content as a string.
        """
        num_wann = self.get_num_wann(structure) if num_wann is None else num_wann
        # NBANDS pre-rounded to a multiple of NPAR (= ntasks) so VASP never
        # changes it, keeping win num_bands identical to the INCAR NBANDS.
        num_bands = self.get_num_bands(num_wann)
        elements = self._projection_elements(structure)

        lines = [
            f"num_wann = {num_wann}",
            f"num_bands = {num_bands}",
            "dis_win_min = -50",
            "dis_win_max = 50",
            "dis_froz_min = -10",
            "dis_froz_max = 10",
            "begin projections",
        ]
        lines.extend(f"{el}: s,p,d" for el in elements)
        lines.extend(
            [
                "end projections",
                "write_hr = .true.",
                "write_xyz = .true.",
                "",
            ]
        )
        return "\n".join(lines)

    def get_input_set(
        self,
        structure: Structure | None = None,
        prev_dir: str | Path | None = None,
        potcar_spec: bool = False,
    ) -> VaspInput:
        """Get the VASP input set for an SCF + Wannier90 calculation.

        Merges the Wannier90-interface INCAR defaults with ``user_incar_settings``
        (user settings win) and writes ``wannier90.win`` as a side effect into the
        current working directory (the job's run directory).

        Args:
            structure: Input structure.
            prev_dir: Previous calculation directory.
            potcar_spec: Use POTCAR.spec symbols instead of POTCAR files.

        Returns:
            VaspInput with the Wannier90 INCAR settings applied.
        """
        if structure is None:
            structure = self.structure
        # Snap the R3-relaxed lattice to the detected spacegroup so wannier90's
        # k-mesh setup sees a clean, symmetric grid (see _symmetrize_structure).
        structure = self._symmetrize_structure(structure)
        num_wann = self.get_num_wann(structure)

        merged_incar = {
            **W90_INCAR_DEFAULTS,
            **self.user_incar_settings,
            "NUM_WANN": num_wann,
            "NBANDS": self.get_num_bands(num_wann),
        }
        original = self.user_incar_settings
        self.user_incar_settings = merged_incar
        try:
            input_set = super().get_input_set(
                structure=structure,
                prev_dir=prev_dir,
                potcar_spec=potcar_spec,
            )
        finally:
            # restore so repeated calls always recompute NUM_WANN / NBANDS
            self.user_incar_settings = original

        # Side effect: pre-write wannier90.win with explicit projections so VASP's
        # Wannier90 interface does not fail with "No projection mode".
        win_string = self.get_wannier90_win_string(structure=structure, num_wann=num_wann)
        Path("wannier90.win").write_text(win_string, encoding="utf-8")
        log.debug(
            f"Wrote wannier90.win (num_wann={num_wann}, "
            f"num_bands={self.get_num_bands(num_wann)}, "
            f"projections={self._projection_elements(structure)})"
        )

        return input_set
