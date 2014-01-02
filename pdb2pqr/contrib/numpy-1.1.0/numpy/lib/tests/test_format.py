r''' Test the .npy file format.

Set up:

    >>> import numpy as np
    >>> from cStringIO import StringIO
    >>> from numpy.lib import format
    >>>
    >>> scalars = [
    ...     np.uint8,
    ...     np.int8,
    ...     np.uint16,
    ...     np.int16,
    ...     np.uint32,
    ...     np.int32,
    ...     np.uint64,
    ...     np.int64,
    ...     np.float32,
    ...     np.float64,
    ...     np.complex64,
    ...     np.complex128,
    ...     object,
    ... ]
    >>>
    >>> basic_arrays = []
    >>>
    >>> for scalar in scalars:
    ...     for endian in '<>':
    ...         dtype = np.dtype(scalar).newbyteorder(endian)
    ...         basic = np.arange(15).astype(dtype)
    ...         basic_arrays.extend([
    ...             np.array([], dtype=dtype),
    ...             np.array(10, dtype=dtype),
    ...             basic,
    ...             basic.reshape((3,5)),
    ...             basic.reshape((3,5)).T,
    ...             basic.reshape((3,5))[::-1,::2],
    ...         ])
    ...
    >>>
    >>> Pdescr = [
    ...     ('x', 'i4', (2,)),
    ...     ('y', 'f8', (2, 2)),
    ...     ('z', 'u1')]
    >>>
    >>>
    >>> PbufferT = [
    ...     ([3,2], [[6.,4.],[6.,4.]], 8),
    ...     ([4,3], [[7.,5.],[7.,5.]], 9),
    ...     ]
    >>>
    >>>
    >>> Ndescr = [
    ...     ('x', 'i4', (2,)),
    ...     ('Info', [
    ...         ('value', 'c16'),
    ...         ('y2', 'f8'),
    ...         ('Info2', [
    ...             ('name', 'S2'),
    ...             ('value', 'c16', (2,)),
    ...             ('y3', 'f8', (2,)),
    ...             ('z3', 'u4', (2,))]),
    ...         ('name', 'S2'),
    ...         ('z2', 'b1')]),
    ...     ('color', 'S2'),
    ...     ('info', [
    ...         ('Name', 'U8'),
    ...         ('Value', 'c16')]),
    ...     ('y', 'f8', (2, 2)),
    ...     ('z', 'u1')]
    >>>
    >>>
    >>> NbufferT = [
    ...     ([3,2], (6j, 6., ('nn', [6j,4j], [6.,4.], [1,2]), 'NN', True), 'cc', ('NN', 6j), [[6.,4.],[6.,4.]], 8),
    ...     ([4,3], (7j, 7., ('oo', [7j,5j], [7.,5.], [2,1]), 'OO', False), 'dd', ('OO', 7j), [[7.,5.],[7.,5.]], 9),
    ...     ]
    >>>
    >>>
    >>> record_arrays = [
    ...     np.array(PbufferT, dtype=np.dtype(Pdescr).newbyteorder('<')),
    ...     np.array(NbufferT, dtype=np.dtype(Ndescr).newbyteorder('<')),
    ...     np.array(PbufferT, dtype=np.dtype(Pdescr).newbyteorder('>')),
    ...     np.array(NbufferT, dtype=np.dtype(Ndescr).newbyteorder('>')),
    ... ]

Test the magic string writing.

    >>> format.magic(1, 0)
    '\x93NUMPY\x01\x00'
    >>> format.magic(0, 0)
    '\x93NUMPY\x00\x00'
    >>> format.magic(255, 255)
    '\x93NUMPY\xff\xff'
    >>> format.magic(2, 5)
    '\x93NUMPY\x02\x05'

Test the magic string reading.

    >>> format.read_magic(StringIO(format.magic(1, 0)))
    (1, 0)
    >>> format.read_magic(StringIO(format.magic(0, 0)))
    (0, 0)
    >>> format.read_magic(StringIO(format.magic(255, 255)))
    (255, 255)
    >>> format.read_magic(StringIO(format.magic(2, 5)))
    (2, 5)

Test the header writing.

    >>> for arr in basic_arrays + record_arrays:
    ...     f = StringIO()
    ...     format.write_array_header_1_0(f, arr)
    ...     print repr(f.getvalue())
    ...
    "F\x00{'descr': '|u1', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '|u1', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '|u1', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '|u1', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '|u1', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '|u1', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '|u1', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '|u1', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '|u1', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '|u1', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '|u1', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '|u1', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '|i1', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '|i1', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '|i1', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '|i1', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '|i1', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '|i1', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '|i1', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '|i1', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '|i1', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '|i1', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '|i1', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '|i1', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '<u2', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '<u2', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '<u2', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '<u2', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '<u2', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '<u2', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '>u2', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '>u2', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '>u2', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '>u2', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '>u2', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '>u2', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '<i2', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '<i2', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '<i2', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '<i2', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '<i2', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '<i2', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '>i2', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '>i2', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '>i2', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '>i2', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '>i2', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '>i2', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '<u4', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '<u4', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '<u4', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '<u4', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '<u4', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '<u4', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '>u4', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '>u4', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '>u4', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '>u4', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '>u4', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '>u4', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '<i4', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '<i4', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '<i4', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '<i4', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '<i4', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '<i4', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '>i4', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '>i4', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '>i4', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '>i4', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '>i4', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '>i4', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '<u8', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '<u8', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '<u8', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '<u8', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '<u8', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '<u8', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '>u8', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '>u8', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '>u8', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '>u8', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '>u8', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '>u8', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '<i8', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '<i8', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '<i8', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '<i8', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '<i8', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '<i8', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '>i8', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '>i8', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '>i8', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '>i8', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '>i8', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '>i8', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '<f4', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '<f4', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '<f4', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '<f4', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '<f4', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '<f4', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '>f4', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '>f4', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '>f4', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '>f4', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '>f4', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '>f4', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '<f8', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '<f8', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '<f8', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '<f8', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '<f8', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '<f8', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '>f8', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '>f8', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '>f8', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '>f8', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '>f8', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '>f8', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '<c8', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '<c8', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '<c8', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '<c8', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '<c8', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '<c8', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '>c8', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '>c8', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '>c8', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '>c8', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '>c8', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '>c8', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '<c16', 'fortran_order': False, 'shape': (0,)}             \n"
    "F\x00{'descr': '<c16', 'fortran_order': False, 'shape': ()}               \n"
    "F\x00{'descr': '<c16', 'fortran_order': False, 'shape': (15,)}            \n"
    "F\x00{'descr': '<c16', 'fortran_order': False, 'shape': (3, 5)}           \n"
    "F\x00{'descr': '<c16', 'fortran_order': True, 'shape': (5, 3)}            \n"
    "F\x00{'descr': '<c16', 'fortran_order': False, 'shape': (3, 3)}           \n"
    "F\x00{'descr': '>c16', 'fortran_order': False, 'shape': (0,)}             \n"
    "F\x00{'descr': '>c16', 'fortran_order': False, 'shape': ()}               \n"
    "F\x00{'descr': '>c16', 'fortran_order': False, 'shape': (15,)}            \n"
    "F\x00{'descr': '>c16', 'fortran_order': False, 'shape': (3, 5)}           \n"
    "F\x00{'descr': '>c16', 'fortran_order': True, 'shape': (5, 3)}            \n"
    "F\x00{'descr': '>c16', 'fortran_order': False, 'shape': (3, 3)}           \n"
    "F\x00{'descr': '|O4', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '|O4', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '|O4', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '|O4', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '|O4', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '|O4', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "F\x00{'descr': '|O4', 'fortran_order': False, 'shape': (0,)}              \n"
    "F\x00{'descr': '|O4', 'fortran_order': False, 'shape': ()}                \n"
    "F\x00{'descr': '|O4', 'fortran_order': False, 'shape': (15,)}             \n"
    "F\x00{'descr': '|O4', 'fortran_order': False, 'shape': (3, 5)}            \n"
    "F\x00{'descr': '|O4', 'fortran_order': True, 'shape': (5, 3)}             \n"
    "F\x00{'descr': '|O4', 'fortran_order': False, 'shape': (3, 3)}            \n"
    "v\x00{'descr': [('x', '<i4', (2,)), ('y', '<f8', (2, 2)), ('z', '|u1')],\n 'fortran_order': False,\n 'shape': (2,)}         \n"
    "\x16\x02{'descr': [('x', '<i4', (2,)),\n           ('Info',\n            [('value', '<c16'),\n             ('y2', '<f8'),\n             ('Info2',\n              [('name', '|S2'),\n               ('value', '<c16', (2,)),\n               ('y3', '<f8', (2,)),\n               ('z3', '<u4', (2,))]),\n             ('name', '|S2'),\n             ('z2', '|b1')]),\n           ('color', '|S2'),\n           ('info', [('Name', '<U8'), ('Value', '<c16')]),\n           ('y', '<f8', (2, 2)),\n           ('z', '|u1')],\n 'fortran_order': False,\n 'shape': (2,)}      \n"
    "v\x00{'descr': [('x', '>i4', (2,)), ('y', '>f8', (2, 2)), ('z', '|u1')],\n 'fortran_order': False,\n 'shape': (2,)}         \n"
    "\x16\x02{'descr': [('x', '>i4', (2,)),\n           ('Info',\n            [('value', '>c16'),\n             ('y2', '>f8'),\n             ('Info2',\n              [('name', '|S2'),\n               ('value', '>c16', (2,)),\n               ('y3', '>f8', (2,)),\n               ('z3', '>u4', (2,))]),\n             ('name', '|S2'),\n             ('z2', '|b1')]),\n           ('color', '|S2'),\n           ('info', [('Name', '>U8'), ('Value', '>c16')]),\n           ('y', '>f8', (2, 2)),\n           ('z', '|u1')],\n 'fortran_order': False,\n 'shape': (2,)}      \n"
'''


from cStringIO import StringIO
import os
import shutil
import tempfile

import numpy as np
from numpy.testing import assert_array_equal, raises

from numpy.lib import format


tempdir = None

# Module-level setup.
def setup_module():
    global tempdir
    tempdir = tempfile.mkdtemp()

def teardown_module():
    global tempdir
    if tempdir is not None and os.path.isdir(tempdir):
        shutil.rmtree(tempdir)
        tempdir = None


# Generate some basic arrays to test with.
scalars = [
    np.uint8,
    np.int8,
    np.uint16,
    np.int16,
    np.uint32,
    np.int32,
    np.uint64,
    np.int64,
    np.float32,
    np.float64,
    np.complex64,
    np.complex128,
    object,
]
basic_arrays = []
for scalar in scalars:
    for endian in '<>':
        dtype = np.dtype(scalar).newbyteorder(endian)
        basic = np.arange(15).astype(dtype)
        basic_arrays.extend([
            # Empty
            np.array([], dtype=dtype),
            # Rank-0
            np.array(10, dtype=dtype),
            # 1-D
            basic,
            # 2-D C-contiguous
            basic.reshape((3,5)),
            # 2-D F-contiguous
            basic.reshape((3,5)).T,
            # 2-D non-contiguous
            basic.reshape((3,5))[::-1,::2],
        ])

# More complicated record arrays.
# This is the structure of the table used for plain objects:
#
# +-+-+-+
# |x|y|z|
# +-+-+-+

# Structure of a plain array description:
Pdescr = [
    ('x', 'i4', (2,)),
    ('y', 'f8', (2, 2)),
    ('z', 'u1')]

# A plain list of tuples with values for testing:
PbufferT = [
    # x     y                  z
    ([3,2], [[6.,4.],[6.,4.]], 8),
    ([4,3], [[7.,5.],[7.,5.]], 9),
    ]


# This is the structure of the table used for nested objects (DON'T PANIC!):
#
# +-+---------------------------------+-----+----------+-+-+
# |x|Info                             |color|info      |y|z|
# | +-----+--+----------------+----+--+     +----+-----+ | |
# | |value|y2|Info2           |name|z2|     |Name|Value| | |
# | |     |  +----+-----+--+--+    |  |     |    |     | | |
# | |     |  |name|value|y3|z3|    |  |     |    |     | | |
# +-+-----+--+----+-----+--+--+----+--+-----+----+-----+-+-+
#

# The corresponding nested array description:
Ndescr = [
    ('x', 'i4', (2,)),
    ('Info', [
        ('value', 'c16'),
        ('y2', 'f8'),
        ('Info2', [
            ('name', 'S2'),
            ('value', 'c16', (2,)),
            ('y3', 'f8', (2,)),
            ('z3', 'u4', (2,))]),
        ('name', 'S2'),
        ('z2', 'b1')]),
    ('color', 'S2'),
    ('info', [
        ('Name', 'U8'),
        ('Value', 'c16')]),
    ('y', 'f8', (2, 2)),
    ('z', 'u1')]

NbufferT = [
    # x     Info                                                color info        y                  z
    #       value y2 Info2                            name z2         Name Value
    #                name   value    y3       z3
    ([3,2], (6j, 6., ('nn', [6j,4j], [6.,4.], [1,2]), 'NN', True), 'cc', ('NN', 6j), [[6.,4.],[6.,4.]], 8),
    ([4,3], (7j, 7., ('oo', [7j,5j], [7.,5.], [2,1]), 'OO', False), 'dd', ('OO', 7j), [[7.,5.],[7.,5.]], 9),
    ]

record_arrays = [
    np.array(PbufferT, dtype=np.dtype(Pdescr).newbyteorder('<')),
    np.array(NbufferT, dtype=np.dtype(Ndescr).newbyteorder('<')),
    np.array(PbufferT, dtype=np.dtype(Pdescr).newbyteorder('>')),
    np.array(NbufferT, dtype=np.dtype(Ndescr).newbyteorder('>')),
]

def roundtrip(arr):
    f = StringIO()
    format.write_array(f, arr)
    f2 = StringIO(f.getvalue())
    arr2 = format.read_array(f2)
    return arr2

def assert_equal(o1, o2):
    assert o1 == o2


def test_roundtrip():
    for arr in basic_arrays + record_arrays:
        arr2 = roundtrip(arr)
        yield assert_array_equal, arr, arr2

def test_memmap_roundtrip():
    for arr in basic_arrays + record_arrays:
        if arr.dtype.hasobject:
            # Skip these since they can't be mmap'ed.
            continue
        # Write it out normally and through mmap.
        nfn = os.path.join(tempdir, 'normal.npy')
        mfn = os.path.join(tempdir, 'memmap.npy')
        fp = open(nfn, 'wb')
        try:
            format.write_array(fp, arr)
        finally:
            fp.close()

        fortran_order = (arr.flags.f_contiguous and not arr.flags.c_contiguous)
        ma = format.open_memmap(mfn, mode='w+', dtype=arr.dtype,
            shape=arr.shape, fortran_order=fortran_order)
        ma[...] = arr
        del ma

        # Check that both of these files' contents are the same.
        fp = open(nfn, 'rb')
        normal_bytes = fp.read()
        fp.close()
        fp = open(mfn, 'rb')
        memmap_bytes = fp.read()
        fp.close()
        yield assert_equal, normal_bytes, memmap_bytes

        # Check that reading the file using memmap works.
        ma = format.open_memmap(nfn, mode='r')
        yield assert_array_equal, ma, arr
        del ma


def test_write_version_1_0():
    f = StringIO()
    arr = np.arange(1)
    # These should pass.
    format.write_array(f, arr, version=(1, 0))
    format.write_array(f, arr)

    # These should all fail.
    bad_versions = [
        (1, 1),
        (0, 0),
        (0, 1),
        (2, 0),
        (2, 2),
        (255, 255),
    ]
    for version in bad_versions:
        try:
            format.write_array(f, arr, version=version)
        except ValueError:
            pass
        else:
            raise AssertionError("we should have raised a ValueError for the bad version %r" % (version,))


bad_version_magic = [
    '\x93NUMPY\x01\x01',
    '\x93NUMPY\x00\x00',
    '\x93NUMPY\x00\x01',
    '\x93NUMPY\x02\x00',
    '\x93NUMPY\x02\x02',
    '\x93NUMPY\xff\xff',
]
malformed_magic = [
    '\x92NUMPY\x01\x00',
    '\x00NUMPY\x01\x00',
    '\x93numpy\x01\x00',
    '\x93MATLB\x01\x00',
    '\x93NUMPY\x01',
    '\x93NUMPY',
    '',
]

def test_read_magic_bad_magic():
    for magic in malformed_magic:
        f = StringIO(magic)
        yield raises(ValueError)(format.read_magic), f

def test_read_version_1_0_bad_magic():
    for magic in bad_version_magic + malformed_magic:
        f = StringIO(magic)
        yield raises(ValueError)(format.read_array), f
