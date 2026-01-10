"""Custom warnings and exceptions for py-sc-fermi."""

import functools
import warnings


class PySCFermiWarning(UserWarning):
	"""Base class for all py-sc-fermi warnings."""
	pass


class DOSOverflowWarning(PySCFermiWarning):
	"""Warning for overflow during DOS carrier concentration calculations.
	
	This typically occurs during iterative Fermi energy solving when
	extreme values cause overflow in exponential calculations.
	"""
	pass


class DefectOverflowWarning(PySCFermiWarning):
	"""Warning for overflow during defect concentration calculations.
	
	This typically occurs during iterative Fermi energy solving when
	extreme values cause overflow in exponential calculations.
	"""
	pass


def catches_numpy_overflow(warning_class):
	"""Decorator that converts numpy overflow RuntimeWarnings to a custom warning.
	
	Args:
		warning_class: The warning class to emit when overflow is detected.
		
	Returns:
		Decorator function.
	"""
	def decorator(func):
		@functools.wraps(func)
		def wrapper(*args, **kwargs):
			with warnings.catch_warnings(record=True) as caught:
				warnings.filterwarnings("always", message="overflow", category=RuntimeWarning)
				result = func(*args, **kwargs)
			
			for w in caught:
				if issubclass(w.category, RuntimeWarning) and "overflow" in str(w.message):
					warnings.warn(
						f"Overflow encountered in {func.__name__}. "
						"This may occur during iterative Fermi energy solving "
						"and can likely be ignored if final results are reasonable.",
						warning_class,
					)
					break  # Only emit once per call
			
			return result
		return wrapper
	return decorator