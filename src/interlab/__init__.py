import copy

from . import metrics
from .project import Project

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _version

try:
    __version__ = _version("interlab_py")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "999"


def fix_spectrum(data):
    """Removes small values from an NMR spectrum and replaces them with an arbitrarily small value, in this case 1.0e-16

    :param anndata: The data to be corrected
    :returns: anndata_fixed, the corrected data
    """
    data_fixed = copy.deepcopy(data)
    eps = 1e-16
    data_fixed[data < eps] = eps
    return data_fixed


__all__ = ["__version__", "metrics", "Project", "fix_spectrum"]
