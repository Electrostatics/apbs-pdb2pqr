"""numerictypes: Define the numeric type objects

This module is designed so 'from numerictypes import *' is safe.
Exported symbols include:

  Dictionary with all registered number types (including aliases):
    typeDict

  Type objects (not all will be available, depends on platform):
      see variable sctypes for which ones you have

    Bit-width names

    int8 int16 int32 int64 int128
    uint8 uint16 uint32 uint64 uint128
    float16 float32 float64 float96 float128 float256
    complex32 complex64 complex128 complex192 complex256 complex512

    c-based names

    bool_

    object_

    void, str_, unicode_

    byte, ubyte,
    short, ushort
    intc, uintc,
    intp, uintp,
    int_, uint,
    longlong, ulonglong,

    single, csingle,
    float_, complex_,
    longfloat, clongfloat,

   As part of the type-hierarchy:    xx -- is bit-width

   generic
     +-> bool_                                  (kind=b)
     +-> number                                 (kind=i)
     |     integer
     |     signedinteger   (intxx)
     |     byte
     |     short
     |     intc
     |     intp           int0
     |     int_
     |     longlong
     +-> unsignedinteger  (uintxx)              (kind=u)
     |     ubyte
     |     ushort
     |     uintc
     |     uintp          uint0
     |     uint_
     |     ulonglong
     +-> inexact
     |   +-> floating           (floatxx)       (kind=f)
     |   |     single
     |   |     float_  (double)
     |   |     longfloat
     |   \-> complexfloating    (complexxx)     (kind=c)
     |         csingle  (singlecomplex)
     |         complex_ (cfloat, cdouble)
     |         clongfloat (longcomplex)
     +-> flexible
     |     character
     |     str_     (string_)                   (kind=S)
     |     unicode_                             (kind=U)
     |     void                                 (kind=V)
     |
     \-> object_ (not used much)                (kind=O)
"""

# we add more at the bottom
__all__ = ['sctypeDict', 'sctypeNA', 'typeDict', 'typeNA', 'sctypes',
           'ScalarType', 'obj2sctype', 'cast', 'nbytes', 'sctype2char',
           'maximum_sctype', 'issctype', 'typecodes', 'find_common_type']

from numpy.core.multiarray import typeinfo, ndarray, array, empty, dtype
import types as _types

# we don't export these for import *, but we do want them accessible
# as numerictypes.bool, etc.
from __builtin__ import bool, int, long, float, complex, object, unicode, str


# String-handling utilities to avoid locale-dependence.

import string
LOWER_TABLE = string.maketrans(string.ascii_uppercase, string.ascii_lowercase)
UPPER_TABLE = string.maketrans(string.ascii_lowercase, string.ascii_uppercase)

def english_lower(s):
    """ Apply English case rules to convert ASCII strings to all lower case.

    This is an internal utility function to replace calls to str.lower() such
    that we can avoid changing behavior with changing locales. In particular,
    Turkish has distinct dotted and dotless variants of the Latin letter "I" in
    both lowercase and uppercase. Thus, "I".lower() != "i" in a "tr" locale.

    Parameters
    ----------
    s : str

    Returns
    -------
    lowered : str

    Examples
    --------
    >>> from numpy.lib.utils import english_lower
    >>> english_lower('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_')
    'abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyz0123456789_'
    >>> english_upper('')
    ''
    """
    lowered = s.translate(LOWER_TABLE)
    return lowered

def english_upper(s):
    """ Apply English case rules to convert ASCII strings to all upper case.

    This is an internal utility function to replace calls to str.upper() such
    that we can avoid changing behavior with changing locales. In particular,
    Turkish has distinct dotted and dotless variants of the Latin letter "I" in
    both lowercase and uppercase. Thus, "i".upper() != "I" in a "tr" locale.

    Parameters
    ----------
    s : str

    Returns
    -------
    uppered : str

    Examples
    --------
    >>> from numpy.lib.utils import english_upper
    >>> english_upper('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_')
    'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'
    >>> english_upper('')
    ''
    """
    uppered = s.translate(UPPER_TABLE)
    return uppered

def english_capitalize(s):
    """ Apply English case rules to convert the first character of an ASCII
    string to upper case.

    This is an internal utility function to replace calls to str.capitalize()
    such that we can avoid changing behavior with changing locales.

    Parameters
    ----------
    s : str

    Returns
    -------
    capitalized : str

    Examples
    --------
    >>> from numpy.lib.utils import english_capitalize
    >>> english_capitalize('int8')
    'Int8'
    >>> english_capitalize('Int8')
    'Int8'
    >>> english_capitalize('')
    ''
    """
    if s:
        return english_upper(s[0]) + s[1:]
    else:
        return s


sctypeDict = {}      # Contains all leaf-node scalar types with aliases
sctypeNA = {}        # Contails all leaf-node types -> numarray type equivalences
allTypes = {}      # Collect the types we will add to the module here

def _evalname(name):
    k = 0
    for ch in name:
        if ch in '0123456789':
            break
        k += 1
    try:
        bits = int(name[k:])
    except ValueError:
        bits = 0
    base = name[:k]
    return base, bits

def bitname(obj):
    """Return a bit-width name for a given type object"""
    name = obj.__name__
    base = ''
    char = ''
    try:
        if name[-1] == '_':
            newname = name[:-1]
        else:
            newname = name
        info = typeinfo[english_upper(newname)]
        assert(info[-1] == obj)  # sanity check
        bits = info[2]

    except KeyError:     # bit-width name
        base, bits = _evalname(name)
        char = base[0]

    if name == 'bool_':
        char = 'b'
        base = 'bool'
    elif name=='string_':
        char = 'S'
        base = 'string'
    elif name=='unicode_':
        char = 'U'
        base = 'unicode'
    elif name=='void':
        char = 'V'
        base = 'void'
    elif name=='object_':
        char = 'O'
        base = 'object'
        bits = 0

    bytes = bits / 8

    if char != '' and bytes != 0:
        char = "%s%d" % (char, bytes)

    return base, bits, char


def _add_types():
    for a in typeinfo.keys():
        name = english_lower(a)
        if isinstance(typeinfo[a], tuple):
            typeobj = typeinfo[a][-1]

            # define C-name and insert typenum and typechar references also
            allTypes[name] = typeobj
            sctypeDict[name] = typeobj
            sctypeDict[typeinfo[a][0]] = typeobj
            sctypeDict[typeinfo[a][1]] = typeobj

        else:  # generic class
            allTypes[name] = typeinfo[a]
_add_types()

def _add_aliases():
    for a in typeinfo.keys():
        name = english_lower(a)
        if not isinstance(typeinfo[a], tuple):
            continue
        typeobj = typeinfo[a][-1]
        # insert bit-width version for this class (if relevant)
        base, bit, char = bitname(typeobj)
        if base[-3:] == 'int' or char[0] in 'ui': continue
        if base != '':
            myname = "%s%d" % (base, bit)
            if (name != 'longdouble' and name != 'clongdouble') or \
                   myname not in allTypes.keys():
                allTypes[myname] = typeobj
                sctypeDict[myname] = typeobj
                if base == 'complex':
                    na_name = '%s%d' % (english_capitalize(base), bit/2)
                elif base == 'bool':
                    na_name = english_capitalize(base)
                    sctypeDict[na_name] = typeobj
                else:
                    na_name = "%s%d" % (english_capitalize(base), bit)
                    sctypeDict[na_name] = typeobj
                sctypeNA[na_name] = typeobj
                sctypeDict[na_name] = typeobj
                sctypeNA[typeobj] = na_name
                sctypeNA[typeinfo[a][0]] = na_name
        if char != '':
            sctypeDict[char] = typeobj
            sctypeNA[char] = na_name
_add_aliases()

# Integers handled so that
# The int32, int64 types should agree exactly with
#  PyArray_INT32, PyArray_INT64 in C
# We need to enforce the same checking as is done
#  in arrayobject.h where the order of getting a
#  bit-width match is:
#       long, longlong, int, short, char
#   for int8, int16, int32, int64, int128

def _add_integer_aliases():
    _ctypes = ['LONG', 'LONGLONG', 'INT', 'SHORT', 'BYTE']
    for ctype in _ctypes:
        val = typeinfo[ctype]
        bits = val[2]
        charname = 'i%d' % (bits/8,)
        ucharname = 'u%d' % (bits/8,)
        intname = 'int%d' % bits
        UIntname = 'UInt%d' % bits
        Intname = 'Int%d' % bits
        uval = typeinfo['U'+ctype]
        typeobj = val[-1]
        utypeobj = uval[-1]
        if intname not in allTypes.keys():
            uintname = 'uint%d' % bits
            allTypes[intname] = typeobj
            allTypes[uintname] = utypeobj
            sctypeDict[intname] = typeobj
            sctypeDict[uintname] = utypeobj
            sctypeDict[Intname] = typeobj
            sctypeDict[UIntname] = utypeobj
            sctypeDict[charname] = typeobj
            sctypeDict[ucharname] = utypeobj
            sctypeNA[Intname] = typeobj
            sctypeNA[UIntname] = utypeobj
            sctypeNA[charname] = typeobj
            sctypeNA[ucharname] = utypeobj
        sctypeNA[typeobj] = Intname
        sctypeNA[utypeobj] = UIntname
        sctypeNA[val[0]] = Intname
        sctypeNA[uval[0]] = UIntname
_add_integer_aliases()

# We use these later
void = allTypes['void']
generic = allTypes['generic']

#
# Rework the Python names (so that float and complex and int are consistent
#                            with Python usage)
#
def _set_up_aliases():
    type_pairs = [('complex_', 'cdouble'),
                  ('int0', 'intp'),
                  ('uint0', 'uintp'),
                  ('single', 'float'),
                  ('csingle', 'cfloat'),
                  ('singlecomplex', 'cfloat'),
                  ('float_', 'double'),
                  ('intc', 'int'),
                  ('uintc', 'uint'),
                  ('int_', 'long'),
                  ('uint', 'ulong'),
                  ('cfloat', 'cdouble'),
                  ('longfloat', 'longdouble'),
                  ('clongfloat', 'clongdouble'),
                  ('longcomplex', 'clongdouble'),
                  ('bool_', 'bool'),
                  ('unicode_', 'unicode'),
                  ('str_', 'string'),
                  ('string_', 'string'),
                  ('object_', 'object')]
    for alias, t in type_pairs:
        allTypes[alias] = allTypes[t]
        sctypeDict[alias] = sctypeDict[t]
    # Remove aliases overriding python types and modules
    for t in ['ulong', 'object', 'unicode', 'int', 'long', 'float',
              'complex', 'bool', 'string']:
        try:
            del allTypes[t]
            del sctypeDict[t]
        except KeyError:
            pass
_set_up_aliases()

# Now, construct dictionary to lookup character codes from types
_sctype2char_dict = {}
def _construct_char_code_lookup():
    for name in typeinfo.keys():
        tup = typeinfo[name]
        if isinstance(tup, tuple):
            if tup[0] not in ['p','P']:
                _sctype2char_dict[tup[-1]] = tup[0]
_construct_char_code_lookup()


sctypes = {'int': [],
           'uint':[],
           'float':[],
           'complex':[],
           'others':[bool,object,str,unicode,void]}

def _add_array_type(typename, bits):
    try:
        t = allTypes['%s%d' % (typename, bits)]
    except KeyError:
        pass
    else:
        sctypes[typename].append(t)

def _set_array_types():
    ibytes = [1, 2, 4, 8, 16, 32, 64]
    fbytes = [2, 4, 8, 10, 12, 16, 32, 64]
    for bytes in ibytes:
        bits = 8*bytes
        _add_array_type('int', bits)
        _add_array_type('uint', bits)
    for bytes in fbytes:
        bits = 8*bytes
        _add_array_type('float', bits)
        _add_array_type('complex', 2*bits)
    _gi = dtype('p')
    if _gi.type not in sctypes['int']:
        indx = 0
        sz = _gi.itemsize
        _lst = sctypes['int']
        while (indx < len(_lst) and sz >= _lst[indx](0).itemsize):
            indx += 1
        sctypes['int'].insert(indx, _gi.type)
        sctypes['uint'].insert(indx, dtype('P').type)
_set_array_types()


genericTypeRank = ['bool', 'int8', 'uint8', 'int16', 'uint16',
                   'int32', 'uint32', 'int64', 'uint64', 'int128',
                   'uint128', 'float16',
                   'float32', 'float64', 'float80', 'float96', 'float128',
                   'float256',
                   'complex32', 'complex64', 'complex128', 'complex160',
                   'complex192', 'complex256', 'complex512', 'object']

def maximum_sctype(t):
    """returns the sctype of highest precision of the same general kind as 't'"""
    g = obj2sctype(t)
    if g is None:
        return t
    t = g
    name = t.__name__
    base, bits = _evalname(name)
    if bits == 0:
        return t
    else:
        return sctypes[base][-1]

_python_types = {int : 'int_',
                 float: 'float_',
                 complex: 'complex_',
                 bool: 'bool_',
                 str: 'string_',
                 unicode: 'unicode_',
                 _types.BufferType: 'void',
                }
def _python_type(t):
    """returns the type corresponding to a certain Python type"""
    if not isinstance(t, _types.TypeType):
        t = type(t)
    return allTypes[_python_types.get(t, 'object_')]

def issctype(rep):
    """Determines whether the given object represents
    a numeric array type."""
    if not isinstance(rep, (type, dtype)):
        return False
    try:
        res = obj2sctype(rep)
        if res and res != object_:
            return True
        return False
    except:
        return False

def obj2sctype(rep, default=None):
    try:
        if issubclass(rep, generic):
            return rep
    except TypeError:
        pass
    if isinstance(rep, dtype):
        return rep.type
    if isinstance(rep, type):
        return _python_type(rep)
    if isinstance(rep, ndarray):
        return rep.dtype.type
    try:
        res = dtype(rep)
    except:
        return default
    return res.type


# This dictionary allows look up based on any alias for an array data-type
class _typedict(dict):
    def __getitem__(self, obj):
        return dict.__getitem__(self, obj2sctype(obj))

nbytes = _typedict()
_alignment = _typedict()
_maxvals = _typedict()
_minvals = _typedict()
def _construct_lookups():
    for name, val in typeinfo.iteritems():
        if not isinstance(val, tuple):
            continue
        obj = val[-1]
        nbytes[obj] = val[2] / 8
        _alignment[obj] = val[3]
        if (len(val) > 5):
            _maxvals[obj] = val[4]
            _minvals[obj] = val[5]
        else:
            _maxvals[obj] = None
            _minvals[obj] = None

_construct_lookups()

def sctype2char(sctype):
    sctype = obj2sctype(sctype)
    if sctype is None:
        raise ValueError, "unrecognized type"
    return _sctype2char_dict[sctype]

# Create dictionary of casting functions that wrap sequences
# indexed by type or type character


cast = _typedict()
ScalarType = [_types.IntType, _types.FloatType,
              _types.ComplexType, _types.LongType, _types.BooleanType,
              _types.StringType, _types.UnicodeType, _types.BufferType]
ScalarType.extend(_sctype2char_dict.keys())
ScalarType = tuple(ScalarType)
for key in _sctype2char_dict.keys():
    cast[key] = lambda x, k=key : array(x, copy=False).astype(k)


_unicodesize = array('u','U1').itemsize

# Create the typestring lookup dictionary
_typestr = _typedict()
for key in _sctype2char_dict.keys():
    if issubclass(key, allTypes['flexible']):
        _typestr[key] = _sctype2char_dict[key]
    else:
        _typestr[key] = empty((1,),key).dtype.str[1:]

# Make sure all typestrings are in sctypeDict
for key, val in _typestr.items():
    if val not in sctypeDict:
        sctypeDict[val] = key

# Add additional strings to the sctypeDict

_toadd = ['int', 'float', 'complex', 'bool', 'object', 'string', ('str', allTypes['string_']),
          'unicode', 'object', ('a', allTypes['string_'])]

for name in _toadd:
    if isinstance(name, tuple):
        sctypeDict[name[0]] = name[1]
    else:
        sctypeDict[name] = allTypes['%s_' % name]

del _toadd, name

# Now add the types we've determined to this module
for key in allTypes:
    globals()[key] = allTypes[key]
    __all__.append(key)

del key

typecodes = {'Character':'c',
             'Integer':'bhilqp',
             'UnsignedInteger':'BHILQP',
             'Float':'fdg',
             'Complex':'FDG',
             'AllInteger':'bBhHiIlLqQpP',
             'AllFloat':'fdgFDG',
             'All':'?bhilqpBHILQPfdgFDGSUVO'}

# backwards compatibility --- deprecated name
typeDict = sctypeDict
typeNA = sctypeNA

# b -> boolean
# u -> unsigned integer
# i -> signed integer
# f -> floating point
# c -> complex
# S -> string
# U -> Unicode string
# V -> record
# O -> Python object
_kind_list = ['b', 'u', 'i', 'f', 'c', 'S', 'U', 'V', 'O']

__test_types = typecodes['AllInteger'][:-2]+typecodes['AllFloat']+'O'
__len_test_types = len(__test_types)

# Keep incrementing until a common type both can be coerced to
#  is found.  Otherwise, return None
def _find_common_coerce(a, b):
    if a > b:
        return a
    try:
        thisind = __test_types.index(a.char)
    except ValueError:
        return None
    while thisind < __len_test_types:
        newdtype = dtype(__test_types[thisind])
        if newdtype >= b and newdtype >= a:
            return newdtype
        thisind += 1
    return None


def find_common_type(array_types, scalar_types):
    """Determine common type following standard coercion rules

    Parameters
    ----------
    array_types : sequence
        A list of dtype convertible objects representing arrays
    scalar_types : sequence
        A list of dtype convertible objects representing scalars

    Returns
    -------
    datatype : dtype
        The common data-type which is the maximum of the array_types
        ignoring the scalar_types unless the maximum of the scalar_types
        is of a different kind.

        If the kinds is not understood, then None is returned.
    """
    array_types = [dtype(x) for x in array_types]
    scalar_types = [dtype(x) for x in scalar_types]

    if len(scalar_types) == 0:
        if len(array_types) == 0:
            return None
        else:
            return max(array_types)
    if len(array_types) == 0:
        return max(scalar_types)

    maxa = max(array_types)
    maxsc = max(scalar_types)

    try:
        index_a = _kind_list.index(maxa.kind)
        index_sc = _kind_list.index(maxsc.kind)
    except ValueError:
        return None

    if index_sc > index_a:
        return _find_common_coerce(maxsc,maxa)
    else:
        return maxa
