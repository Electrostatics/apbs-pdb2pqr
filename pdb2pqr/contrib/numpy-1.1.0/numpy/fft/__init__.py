# To get sub-modules
from info import __doc__

from fftpack import *
from helper import *

def test(level=1, verbosity=1):
    from numpy.testing import NumpyTest
    return NumpyTest().test(level, verbosity)
