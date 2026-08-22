"""
HT-VASP - Endmember model
"""

import logging
import os
from pathlib import Path

from dotenv import load_dotenv
from pymatgen.core import Element, Structure

from htvasp.model.templates import TemplateRegistry

log = logging.getLogger(__name__)

load_dotenv(Path(__file__).resolve().parents[2] / ".env")
api_key = os.getenv("MP_API_KEY")

class Endmember:
    def __init__(
            self,
            api_key: str | None = api_key,
            phase_func: dict| None = None,
            elem_mpid: dict | None = None,
            registry: TemplateRegistry | None = None,
        ):
        self.api_key = api_key
        self.registry = registry or TemplateRegistry()
        self.phase_func = phase_func or {
            "SER": self.get_ser,
            "BCC": self.get_bcc,
            "FCC": self.get_fcc,
            "HCP": self.get_hcp,
        }
        self.elem_mpid = elem_mpid or {
            "Al": "mp-134",     # FCC
            "Co": "mp-54",      # HCP
            "Cr": "mp-90",      # BCC
            "Cu": "mp-30",      # FCC
            "Fe": "mp-13",      # BCC
            "Mn": "mp-35",      # CBCC
            "Nb": "mp-75",      # BCC
            "Ni": "mp-23",      # FCC
            "Ta": "mp-50",      # BCC
            "Ti": "mp-46",      # HCP
            "V": "mp-146",      # BCC
            "W": "mp-91",       # BCC
            "Zr": "mp-131",     # HCP
        }

    def get_structure(self, mpid: str) -> Structure:
        """Fetches the structure from the Materials Project with mpid.

        Args:
            mpid: Materials Project ID
        Returns:
            Structure
        """
        if not self.api_key:
            raise ValueError("MP_API_KEY is required for Materials Project structure lookup")
        try:
            from mp_api.client import MPRester
        except ImportError as exc:
            raise ImportError("mp-api is required for Materials Project structure lookup") from exc
        with MPRester(self.api_key) as mpr:
            structure = mpr.get_structure_by_material_id(mpid)
        return structure[0] if isinstance(structure, list) else structure

    def get_ser(self, element: str) -> Structure:
        """
        Fetches the structure SER of element and mpid.

        Args:
            element: Element to fetch structure for
        Returns:
            Structure SER
        """
        mpid = self.elem_mpid[element]
        struct = self.get_structure(mpid)
        log.info(f"Fetched SER for {element}: {mpid}")
        return struct

    def get_bcc(self, elem1: str, elem2: str) -> Structure:
        """
        Creates a BCC alloy structure from two elements.
        Template is Ni_1a-Al_1b, Pm-3m, mp-1487.

        Args:
            elem1: First element symbol
            elem2: Second element symbol
        Returns:
            Structure
        """
        templ = self.registry.get("BCC")
        source = self.registry.get_record("BCC").source_species
        struct = templ.replace_species({source[0]: Element(elem1), source[1]: Element(elem2)})
        return struct.get_sorted_structure(lambda site: str(site.specie).lower())

    def get_fcc(self, elem1: str, elem2: str) -> Structure:
        """
        Creates a FCC alloy structure from two elements.
        Template is Au_1a-Cu_3c, Pm-3m, mp-2258.

        Args:
            elem1: First element symbol
            elem2: Second element symbol
        Returns:
            Structure
        """
        templ = self.registry.get("FCC")
        source = self.registry.get_record("FCC").source_species
        struct = templ.replace_species({source[0]: Element(elem1), source[1]: Element(elem2)})
        return struct.get_sorted_structure(lambda site: str(site.specie).lower())

    def get_hcp(self, elem1: str, elem2: str) -> Structure:
        """
        Creates a HCP alloy structure from two elements.
        Template is Sn_2c-Ni_6h, P63/mmc, mp-20112.

        Args:
            elem1: First element symbol
            elem2: Second element symbol
        Returns:
            Structure
        """
        templ = self.registry.get("HCP")
        source = self.registry.get_record("HCP").source_species
        struct = templ.replace_species({source[0]: Element(elem1), source[1]: Element(elem2)})
        return struct.get_sorted_structure(lambda site: str(site.specie).lower())

    def get_poscar(self, name: str, outdir: Path | str) -> Structure:
        """
        Reads a POSCAR file and returns the structure.

        Args:
            name: Name of the POSCAR file
            outdir: Directory where the POSCAR file is located
        Returns:
            Structure
        """
        outdir = Path(outdir) if isinstance(outdir, str) else outdir
        outdir.mkdir(parents=True, exist_ok=True)
        poscar = outdir.joinpath(f"{name}.vasp")
        if poscar.is_file():
            return Structure.from_file(poscar)
        else:
            parts = name.split("-")
            phase = parts[0]
            elems = parts[1:]
            struct = self.phase_func[phase](*elems)
            struct.to(fmt="poscar", filename=poscar)
            return struct
