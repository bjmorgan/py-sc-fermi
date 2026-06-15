"""Warning categories and numpy overflow handling for py-sc-fermi."""

from collections.abc import Callable
from functools import wraps
from typing import ParamSpec, TypeVar

import numpy as np


class DiluteLimitWarning(UserWarning):
    """Warns that a defect's solved site occupancy is high enough that
    py-sc-fermi's dilute, non-interacting-defect assumption may no longer hold,
    so the results may be non-physical.

    Emitted by :meth:`DefectSystem.get_sc_fermi` -- and therefore by every
    results path that solves (``report``, ``concentration_dict``,
    ``site_percentages``) -- when any species' occupancy exceeds
    ``DefectSystem.occupancy_warning_threshold``, at most once per
    ``DefectSystem`` instance. A dedicated ``UserWarning`` subclass so it can be
    filtered independently of other warnings.
    """


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