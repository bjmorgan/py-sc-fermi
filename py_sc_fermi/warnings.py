"""Warning categories and numpy overflow handling for py-sc-fermi."""

from collections.abc import Callable
from functools import wraps
from typing import ParamSpec, TypeVar

import numpy as np


class PyScFermiWarning(UserWarning):
    """Base class for py-sc-fermi warnings, so callers can filter them as a group."""


class DiluteLimitWarning(PyScFermiWarning):
    """Warns that a defect's solved site occupancy is high enough that
    py-sc-fermi's dilute, non-interacting-defect assumption may no longer hold,
    so the results may be non-physical.

    Emitted by :meth:`DefectSystem.get_sc_fermi` -- and therefore by every
    results path that solves (``result``, ``site_percentages``,
    ``element_chemical_potential_shifts``) -- when any species' occupancy
    exceeds ``DefectSystem.occupancy_warning_threshold``, at most once per
    ``DefectSystem`` instance. A dedicated ``PyScFermiWarning`` subclass, so it
    can be filtered on its own or together with other py-sc-fermi warnings.
    """


class UnrecognisedKeyWarning(PyScFermiWarning):
    """Warns that a ``from_dict`` input contained keys that were ignored."""


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