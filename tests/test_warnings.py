import subprocess
import sys
import unittest
import warnings


class TestWarningsNotSuppressed(unittest.TestCase):
	"""Test that importing py_sc_fermi does not suppress unrelated warnings."""

	def test_user_warning_not_silenced(self):
		"""UserWarning should not be suppressed after importing py_sc_fermi."""
		script = """
import warnings
from py_sc_fermi.defect_system import DefectSystem
warnings.warn("test user warning", UserWarning)
"""
		result = subprocess.run(
			[sys.executable, "-c", script],
			capture_output=True,
			text=True,
		)
		self.assertIn("UserWarning", result.stderr)
		self.assertIn("test user warning", result.stderr)

	def test_deprecation_warning_not_silenced(self):
			"""DeprecationWarning should not be suppressed after importing py_sc_fermi."""
			script = """
import warnings
warnings.filterwarnings("always", category=DeprecationWarning)
from py_sc_fermi.defect_system import DefectSystem
warnings.warn("test deprecation warning", DeprecationWarning)
"""
			result = subprocess.run(
				[sys.executable, "-c", script],
				capture_output=True,
				text=True,
			)
			self.assertIn("DeprecationWarning", result.stderr)
			self.assertIn("test deprecation warning", result.stderr)
	
	def test_runtime_warning_not_intercepted(self):
		"""RuntimeWarning should go through normal warning machinery."""
		script = """
import warnings
from py_sc_fermi.defect_system import DefectSystem
warnings.warn("test runtime warning", RuntimeWarning)
"""
		result = subprocess.run(
			[sys.executable, "-c", script],
			capture_output=True,
			text=True,
		)
		self.assertIn("RuntimeWarning", result.stderr)
		self.assertIn("test runtime warning", result.stderr)
			
		
class TestWarningClasses(unittest.TestCase):
	"""Test custom warning class hierarchy."""

	def test_warning_hierarchy(self):
		"""Custom warnings should inherit from PySCFermiWarning."""
		from py_sc_fermi.warnings import (
			PySCFermiWarning,
			DOSOverflowWarning,
			DefectOverflowWarning,
		)
		
		self.assertTrue(issubclass(PySCFermiWarning, UserWarning))
		self.assertTrue(issubclass(DOSOverflowWarning, PySCFermiWarning))
		self.assertTrue(issubclass(DefectOverflowWarning, PySCFermiWarning))


class TestCatchesNumpyOverflowDecorator(unittest.TestCase):
	"""Test the catches_numpy_overflow decorator."""

	def test_emits_custom_warning_on_overflow(self):
		"""Decorator should emit custom warning when numpy overflow occurs."""
		import numpy as np
		from py_sc_fermi.warnings import catches_numpy_overflow, DOSOverflowWarning

		@catches_numpy_overflow(DOSOverflowWarning)
		def cause_overflow():
			return np.exp(1000)

		with self.assertWarns(DOSOverflowWarning):
			cause_overflow()

	def test_returns_function_result(self):
		"""Decorator should return the wrapped function's result."""
		from py_sc_fermi.warnings import catches_numpy_overflow, DOSOverflowWarning

		@catches_numpy_overflow(DOSOverflowWarning)
		def add(a, b):
			return a + b

		self.assertEqual(add(2, 3), 5)

	def test_no_warning_when_no_overflow(self):
		"""Decorator should not emit warning when no overflow occurs."""
		import numpy as np
		from py_sc_fermi.warnings import catches_numpy_overflow, DOSOverflowWarning

		@catches_numpy_overflow(DOSOverflowWarning)
		def safe_exp():
			return np.exp(1)

		with warnings.catch_warnings(record=True) as caught:
			warnings.simplefilter("always")
			result = safe_exp()

		overflow_warnings = [w for w in caught if issubclass(w.category, DOSOverflowWarning)]
		self.assertEqual(len(overflow_warnings), 0)
		self.assertAlmostEqual(result, np.e)	

if __name__ == "__main__":
	unittest.main()