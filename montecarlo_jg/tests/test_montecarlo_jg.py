"""
Unit and regression test for the montecarlo_jg package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import montecarlo_jg


def test_montecarlo_jg_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "montecarlo_jg" in sys.modules
