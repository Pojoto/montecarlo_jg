
# Import package, test suite, and other packages as needed
import sys

import pytest

import montecarlo_jg


def test_canvas():
    quote = montecarlo_jg.canvas(False)
    assert("The code is but a canvas to our imagination." == quote)
    quote = montecarlo_jg.canvas()
    assert("The code is but a canvas to our imagination.\n\t- Adapted from Henry David Thoreau" == quote)


def test_zen():
    quote = montecarlo_jg.zen(False)
    assert(len(quote) == 894)
    quote = montecarlo_jg.zen(True)
    assert(len(quote) == 906)