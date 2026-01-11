"""Numpy overflow handling for py-sc-fermi."""

from functools import wraps

import numpy as np


def suppresses_numpy_overflow(func):
    """Decorator that suppresses numpy overflow during calculation.

    Overflow commonly occurs during Fermi energy solving when evaluating
    carrier and defect concentrations at extreme energies. This is expected
    behaviour and the results (inf) are mathematically correct.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        old_settings = np.seterr(over='ignore')
        try:
            return func(*args, **kwargs)
        finally:
            np.seterr(**old_settings)
    return wrapper