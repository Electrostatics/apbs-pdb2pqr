# Functions that already have the correct syntax or miscellaneous functions


__all__ = ['sort', 'copy_reg', 'clip', 'rank',
           'sign', 'shape', 'types', 'allclose', 'size',
           'choose', 'swapaxes', 'array_str',
           'pi', 'math', 'concatenate', 'putmask', 'put',
           'around', 'vdot', 'transpose', 'array2string', 'diagonal',
           'searchsorted', 'copy', 'resize',
           'array_repr', 'e', 'StringIO', 'pickle',
           'argsort', 'convolve', 'cross_correlate',
           'dot', 'outerproduct', 'innerproduct', 'insert']

import types
import StringIO
import pickle
import math
import copy
import copy_reg

from numpy import sort, clip, rank, sign, shape, putmask, allclose, size,\
     choose, swapaxes, array_str, array_repr, e, pi, put, \
     resize, around, concatenate, vdot, transpose, \
     diagonal, searchsorted, argsort, convolve, dot, \
     outer as outerproduct, inner as innerproduct, \
     correlate as cross_correlate, \
     place as insert

from array_printer import array2string
