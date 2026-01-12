import subprocess
import sys
import unittest
import warnings

import numpy as np


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

class TestWarnings(unittest.TestCase):
        
    def test_suppresses_numpy_overflow_decorator(self):
        """Decorator should suppress numpy overflow warnings."""
        from py_sc_fermi.warnings import suppresses_numpy_overflow
        
        @suppresses_numpy_overflow
        def cause_overflow():
            return np.exp(1000)
        
        with warnings.catch_warnings(record=True) as caught:
            warnings.filterwarnings("always")
            result = cause_overflow()
        
        assert result == np.inf
        assert len(caught) == 0

if __name__ == "__main__":
    unittest.main()