"""
HT-VASP - Endmember model
"""

import logging
import os
from pathlib import Path

from dotenv import load_dotenv
from mp_api.client import MPRester
from pymatgen.core import Element, Structure

log = logging.getLogger(__name__)

load_dotenv()

api_key = os.getenv("MP_API_KEY", "")
if not api_key:
    raise ValueError("MP_API_KEY not set in .env file or environment variables")

class Endmember:
    def __init__(
            self, 
            api_key: str = api_key, 
            phase_func: dict| None = None,
            elem_mpid: dict | None = None,
        ):
        self.api_key = api_key
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
        templ = self.get_structure("mp-1487")
        struct = templ.replace_species({"Ni": Element(elem1), "Al": Element(elem2)})
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
        templ = self.get_structure("mp-2258")
        struct = templ.replace_species({"Au": Element(elem1), "Cu": Element(elem2)})
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
        templ = self.get_structure("mp-20112")
        struct = templ.replace_species({"Sn": Element(elem1), "Ni": Element(elem2)})
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
