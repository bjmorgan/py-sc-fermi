"""Numpy overflow handling for py-sc-fermi."""

from collections.abc import Callable
from functools import wraps
from typing import ParamSpec, TypeVar

import numpy as np

P = ParamSpec("P")
R = TypeVar("R")


def suppresses_numpy_overflow(func: Callable[P, R]) -> Callable[P, R]:
    """Decorator that suppresses numpy overflow during calculation.

    Overflow commonly occurs during Fermi energy solving when evaluating
    carrier and defect concentrations at extreme energies. This is expected
    behaviour and the results (inf) are mathematically correct.
    """
    @wraps(func)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        old_settings = np.seterr(over='ignore')
        try:
            return func(*args, **kwargs)
        finally:
            np.seterr(**old_settings)
    return wrapper