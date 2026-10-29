"""Minimal netCDF4 shim so TB2J/HamiltonIO imports without the real package.

Only the EPW code path touches netCDF4 and TB2J's ``wann2J`` never calls it.
The real package is absent from the ht-vasp conda env, so
``python -m TB2J.scripts.wann2J`` died at import time.
"""
import importlib.machinery as _m
import sys
import types

if "netCDF4" not in sys.modules:
    mod = types.ModuleType("netCDF4")

    class Dataset:  # noqa: D401 - inert placeholder
        def __init__(self, *a, **k):
            raise RuntimeError("netCDF4 is not installed; the EPW path is unavailable")

    mod.Dataset = Dataset
    mod.__version__ = "0.0-stub"
    mod.__spec__ = _m.ModuleSpec("netCDF4", loader=None)
    sys.modules["netCDF4"] = mod
