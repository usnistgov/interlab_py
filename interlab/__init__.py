import copy
import sys
import os

sys.path.append(os.path.dirname(__file__))

import metrics
from project import Project

def fix_spectrum(data):
    """Removes small values from an NMR spectrum and replaces them with an arbitrarily small value, in this case 1.0e-16
    
    :param anndata: The data to be corrected
    :returns: anndata_fixed, the corrected data
    """
    data_fixed = copy.deepcopy(data)
    eps = 1e-16
    data_fixed[data < eps] = eps
    return data_fixed

