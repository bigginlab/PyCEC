"""
PyCEC: A ... package for Python
================================================

Documentation is available in the docstrings and
online at <link>

Subpackages
-----------
Using any of these subpackages requires an explicit import. For example,
``import scipy.cluster``.

::

 cluster                      --- Vector Quantization / Kmeans
 constants                    --- Physical and mathematical constants and units
 datasets                     --- Dataset methods


Public API in the main PyCEC namespace
--------------------------------------
::

 __version__       --- PyCEC version string
 LowLevelCallable  --- Low-level callback function
 show_config       --- Show PyCEC build configuration
 test              --- Run PyCEC unittests
"""

#export PYTHONPATH="${PYTHONPATH}:/biggin/b222/catz0163/pept/dynamics/pept_holo/pept_AF_H87P_D342P_v2/qmmm/analysis/PyCEC/PyCEC"

# # For relative imports to work in Python 3.6
# import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from pathlib import Path

def _get_version() -> str:
    """Read VERSION.txt and return its contents."""
    path = Path(__file__).parent.resolve()
    version_file = path / "VERSION.txt"
    return version_file.read_text(encoding="utf-8").strip()


__version__ = _get_version()