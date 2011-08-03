from info import __doc__
from numpy.version import version as __version__

from type_check import *
from index_tricks import *
from function_base import *
from shape_base import *
from twodim_base import *
from ufunclike import *

import scimath as emath
from polynomial import *
from machar import *
from getlimits import *
#import convertcode
from utils import *
from arraysetops import *
from io import *
from financial import *
import math

__all__ = ['emath','math']
__all__ += type_check.__all__
__all__ += index_tricks.__all__
__all__ += function_base.__all__
__all__ += shape_base.__all__
__all__ += twodim_base.__all__
__all__ += ufunclike.__all__
__all__ += polynomial.__all__
__all__ += machar.__all__
__all__ += getlimits.__all__
__all__ += utils.__all__
__all__ += arraysetops.__all__
__all__ += io.__all__
__all__ += financial.__all__

def test(level=1, verbosity=1):
    from numpy.testing import NumpyTest
    return NumpyTest().test(level, verbosity)
