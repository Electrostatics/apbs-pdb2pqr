import re

Zero = "PyUFunc_Zero"
One = "PyUFunc_One"
None_ = "PyUFunc_None"

class TypeDescription(object):
    """Type signature for a ufunc

    Attributes
    ----------

    type: character representing the type
    func_data:
    in_:
    out:
    """
    def __init__(self, type, f=None, in_=None, out=None):
        self.type = type
        self.func_data = f
        if in_ is not None:
            in_ = in_.replace('.', type)
        self.in_ = in_
        if out is not None:
            out = out.replace('.', type)
        self.out = out

    def finish_signature(self, nin, nout):
        if self.in_ is None:
            self.in_ = self.type * nin
        assert len(self.in_) == nin
        if self.out is None:
            self.out = self.type * nout
        assert len(self.out) == nout

_fdata_map = dict(f='%sf', d='%s', g='%sl',
                  F='nc_%sf', D='nc_%s', G='nc_%sl')
def build_func_data(types, f):
    func_data = []
    for t in types:
        d = _fdata_map.get(t, '%s') % (f,)
        func_data.append(d)
    return func_data

def TD(types, f=None, in_=None, out=None):
    if f is not None:
        if isinstance(f, str):
            func_data = build_func_data(types, f)
        else:
            assert len(f) == len(types)
            func_data = f
    else:
        func_data = (None,) * len(types)
    if isinstance(in_, str):
        in_ = (in_,) * len(types)
    elif in_ is None:
        in_ = (None,) * len(types)
    if isinstance(out, str):
        out = (out,) * len(types)
    elif out is None:
        out = (None,) * len(types)
    tds = []
    for t, fd, i, o in zip(types, func_data, in_, out):
        tds.append(TypeDescription(t, f=fd, in_=i, out=o))
    return tds

class Ufunc(object):
    """Description of a ufunc.

    Attributes
    ----------

    nin: number of input arguments
    nout: number of output arguments
    identity: identity element for a two-argument function
    docstring: docstring for the ufunc
    type_descriptions: list of TypeDescription objects
    """
    def __init__(self, nin, nout, identity, docstring,
                 *type_descriptions):
        self.nin = nin
        self.nout = nout
        if identity is None:
            identity = None_
        self.identity = identity
        self.docstring = docstring
        self.type_descriptions = []
        for td in type_descriptions:
            self.type_descriptions.extend(td)
        for td in self.type_descriptions:
            td.finish_signature(self.nin, self.nout)

#each entry in defdict is a Ufunc object.

#name: [string of chars for which it is defined,
#       string of characters using func interface,
#       tuple of strings giving funcs for data,
#       (in, out), or (instr, outstr) giving the signature as character codes,
#       identity,
#       docstring,
#       output specification (optional)
#       ]

chartoname = {'?': 'bool',
              'b': 'byte',
              'B': 'ubyte',
              'h': 'short',
              'H': 'ushort',
              'i': 'int',
              'I': 'uint',
              'l': 'long',
              'L': 'ulong',
              'q': 'longlong',
              'Q': 'ulonglong',
              'f': 'float',
              'd': 'double',
              'g': 'longdouble',
              'F': 'cfloat',
              'D': 'cdouble',
              'G': 'clongdouble',
              'O': 'OBJECT',
              # M is like O, but calls a method of the object instead
              # of a function
              'M': 'OBJECT',
              }

all = '?bBhHiIlLqQfdgFDGO'
O = 'O'
M = 'M'
ints = 'bBhHiIlLqQ'
intsO = ints + O
bints = '?' + ints
bintsO = bints + O
flts = 'fdg'
fltsO = flts + O
fltsM = flts + M
cmplx = 'FDG'
cmplxO = cmplx + O
cmplxM = cmplx + M
inexact = flts + cmplx
noint = inexact+O
nointM = inexact+M
allM = bints+flts+cmplxM
nobool = all[1:]
nobool_or_obj = all[1:-1]
intflt = ints+flts
intfltcmplx = nobool_or_obj
nocmplx = bints+flts
nocmplxO = nocmplx+O
nocmplxM = nocmplx+M
noobj = all[:-1]

defdict = {
'add' :
    Ufunc(2, 1, Zero,
          'adds the arguments elementwise.',
          TD(noobj),
          TD(O, f='PyNumber_Add'),
          ),
'subtract' :
    Ufunc(2, 1, Zero,
          'subtracts the arguments elementwise.',
          TD(noobj),
          TD(O, f='PyNumber_Subtract'),
          ),
'multiply' :
    Ufunc(2, 1, One,
          'multiplies the arguments elementwise.',
          TD(noobj),
          TD(O, f='PyNumber_Multiply'),
          ),
'divide' :
    Ufunc(2, 1, One,
          'divides the arguments elementwise.',
          TD(intfltcmplx),
          TD(O, f='PyNumber_Divide'),
          ),
'floor_divide' :
    Ufunc(2, 1, One,
          'floor divides the arguments elementwise.',
          TD(intfltcmplx),
          TD(O, f='PyNumber_FloorDivide'),
          ),
'true_divide' :
    Ufunc(2, 1, One,
          'true divides the arguments elementwise.',
          TD('bBhH', out='f'),
          TD('iIlLqQ', out='d'),
          TD(flts+cmplx),
          TD(O, f='PyNumber_TrueDivide'),
          ),
'conjugate' :
    Ufunc(1, 1, None,
          'takes the conjugate of x elementwise.',
          TD(nobool_or_obj),
          TD(M, f='conjugate'),
          ),
'fmod' :
    Ufunc(2, 1, Zero,
          'computes (C-like) x1 % x2 elementwise.',
          TD(ints),
          TD(flts, f='fmod'),
          TD(M, f='fmod'),
          ),
'square' :
    Ufunc(1, 1, None,
          'compute x**2.',
          TD(nobool_or_obj),
          TD(O, f='Py_square'),
          ),
'reciprocal' :
    Ufunc(1, 1, None,
          'compute 1/x',
          TD(nobool_or_obj),
          TD(O, f='Py_reciprocal'),
          ),
'ones_like' :
    Ufunc(1, 1, None,
          'returns an array of ones of the shape and typecode of x.',
          TD(nobool_or_obj),
          TD(O, f='Py_get_one'),
          ),
'power' :
    Ufunc(2, 1, One,
          'computes x1**x2 elementwise.',
          TD(ints),
          TD(inexact, f='pow'),
          TD(O, f='PyNumber_Power'),
          ),
'absolute' :
    Ufunc(1, 1, None,
          'takes |x| elementwise.',
          TD(nocmplx),
          TD(cmplx, out=('f', 'd', 'g')),
          TD(O, f='PyNumber_Absolute'),
          ),
'negative' :
    Ufunc(1, 1, None,
          'determines -x elementwise',
          TD(nocmplx),
          TD(cmplx, f='neg'),
          TD(O, f='PyNumber_Negative'),
          ),
'sign' :
    Ufunc(1, 1, None,
          'returns -1 if x < 0 and 0 if x==0 and 1 if x > 0',
          TD(nobool),
          ),
'greater' :
    Ufunc(2, 1, None,
          'returns elementwise x1 > x2 in a bool array.',
          TD(all, out='?'),
          ),
'greater_equal' :
    Ufunc(2, 1, None,
          'returns elementwise x1 >= x2 in a bool array.',
          TD(all, out='?'),
          ),
'less' :
    Ufunc(2, 1, None,
          'returns elementwise x1 < x2 in a bool array.',
          TD(all, out='?'),
          ),
'less_equal' :
    Ufunc(2, 1, None,
          'returns elementwise x1 <= x2 in a bool array',
          TD(all, out='?'),
          ),
'equal' :
    Ufunc(2, 1, None,
          'returns elementwise x1 == x2 in a bool array',
          TD(all, out='?'),
          ),
'not_equal' :
    Ufunc(2, 1, None,
          'returns elementwise x1 |= x2',
          TD(all, out='?'),
          ),
'logical_and' :
    Ufunc(2, 1, One,
          'returns x1 and x2 elementwise.',
          TD(noobj, out='?'),
          TD(M, f='logical_and'),
          ),
'logical_not' :
    Ufunc(1, 1, None,
          'returns not x elementwise.',
          TD(noobj, out='?'),
          TD(M, f='logical_not'),
          ),
'logical_or' :
    Ufunc(2, 1, Zero,
          'returns x1 or x2 elementwise.',
          TD(noobj, out='?'),
          TD(M, f='logical_or'),
          ),
'logical_xor' :
    Ufunc(2, 1, None,
          'returns x1 xor x2 elementwise.',
          TD(noobj, out='?'),
          TD(M, f='logical_xor'),
          ),
'maximum' :
    Ufunc(2, 1, None,
          'returns maximum (if x1 > x2: x1;  else: x2) elementwise.',
          TD(noobj),
          TD(O, f='_npy_ObjectMax')
          ),
'minimum' :
    Ufunc(2, 1, None,
          'returns minimum (if x1 < x2: x1;  else: x2) elementwise',
          TD(noobj),
          TD(O, f='_npy_ObjectMin')
          ),
'bitwise_and' :
    Ufunc(2, 1, One,
          'computes x1 & x2 elementwise.',
          TD(bints),
          TD(O, f='PyNumber_And'),
          ),
'bitwise_or' :
    Ufunc(2, 1, Zero,
          'computes x1 | x2 elementwise.',
          TD(bints),
          TD(O, f='PyNumber_Or'),
          ),
'bitwise_xor' :
    Ufunc(2, 1, None,
          'computes x1 ^ x2 elementwise.',
          TD(bints),
          TD(O, f='PyNumber_Xor'),
          ),
'invert' :
    Ufunc(1, 1, None,
          'computes ~x (bit inversion) elementwise.',
          TD(bints),
          TD(O, f='PyNumber_Invert'),
          ),
'left_shift' :
    Ufunc(2, 1, None,
          'computes x1 << x2 (x1 shifted to left by x2 bits) elementwise.',
          TD(ints),
          TD(O, f='PyNumber_Lshift'),
          ),
'right_shift' :
    Ufunc(2, 1, None,
          'computes x1 >> x2 (x1 shifted to right by x2 bits) elementwise.',
          TD(ints),
          TD(O, f='PyNumber_Rshift'),
          ),
'degrees' :
    Ufunc(1, 1, None,
          'converts angle from radians to degrees',
          TD(fltsM, f='degrees'),
          ),
'radians' :
    Ufunc(1, 1, None,
          'converts angle from degrees to radians',
          TD(fltsM, f='radians'),
          ),
'arccos' :
    Ufunc(1, 1, None,
          'inverse cosine elementwise.',
          TD(inexact, f='acos'),
          TD(M, f='arccos'),
          ),
'arccosh' :
    Ufunc(1, 1, None,
          'inverse hyperbolic cosine elementwise.',
          TD(inexact, f='acosh'),
          TD(M, f='arccosh'),
          ),
'arcsin' :
    Ufunc(1, 1, None,
          'inverse sine elementwise.',
          TD(inexact, f='asin'),
          TD(M, f='arcsin'),
          ),
'arcsinh' :
    Ufunc(1, 1, None,
          'inverse hyperbolic sine elementwise.',
          TD(inexact, f='asinh'),
          TD(M, f='arcsinh'),
          ),
'arctan' :
    Ufunc(1, 1, None,
          'inverse tangent elementwise.',
          TD(inexact, f='atan'),
          TD(M, f='arctan'),
          ),
'arctanh' :
    Ufunc(1, 1, None,
          'inverse hyperbolic tangent elementwise.',
          TD(inexact, f='atanh'),
          TD(M, f='arctanh'),
          ),
'cos' :
    Ufunc(1, 1, None,
          'cosine elementwise.',
          TD(inexact, f='cos'),
          TD(M, f='cos'),
          ),
'sin' :
    Ufunc(1, 1, None,
          'sine elementwise.',
          TD(inexact, f='sin'),
          TD(M, f='sin'),
          ),
'tan' :
    Ufunc(1, 1, None,
          'tangent elementwise.',
          TD(inexact, f='tan'),
          TD(M, f='tan'),
          ),
'cosh' :
    Ufunc(1, 1, None,
          'hyperbolic cosine elementwise.',
          TD(inexact, f='cosh'),
          TD(M, f='cosh'),
          ),
'sinh' :
    Ufunc(1, 1, None,
          'hyperbolic sine elementwise.',
          TD(inexact, f='sinh'),
          TD(M, f='sinh'),
          ),
'tanh' :
    Ufunc(1, 1, None,
          'hyperbolic tangent elementwise.',
          TD(inexact, f='tanh'),
          TD(M, f='tanh'),
          ),
'exp' :
    Ufunc(1, 1, None,
          'e**x elementwise.',
          TD(inexact, f='exp'),
          TD(M, f='exp'),
          ),
'expm1' :
    Ufunc(1, 1, None,
          'e**x-1 elementwise.',
          TD(inexact, f='expm1'),
          TD(M, f='expm1'),
          ),
'log' :
    Ufunc(1, 1, None,
          'logarithm base e elementwise.',
          TD(inexact, f='log'),
          TD(M, f='log'),
          ),
'log10' :
    Ufunc(1, 1, None,
          'logarithm base 10 elementwise.',
          TD(inexact, f='log10'),
          TD(M, f='log10'),
          ),
'log1p' :
    Ufunc(1, 1, None,
          'log(1+x) to base e elementwise.',
          TD(inexact, f='log1p'),
          TD(M, f='log1p'),
          ),
'sqrt' :
    Ufunc(1, 1, None,
          'square-root elementwise. For real x, the domain is restricted to x>=0.',
          TD(inexact, f='sqrt'),
          TD(M, f='sqrt'),
          ),
'ceil' :
    Ufunc(1, 1, None,
          'elementwise smallest integer >= x.',
          TD(flts, f='ceil'),
          TD(M, f='ceil'),
          ),
'fabs' :
    Ufunc(1, 1, None,
          'absolute values.',
          TD(flts, f='fabs'),
          TD(M, f='fabs'),
       ),
'floor' :
    Ufunc(1, 1, None,
          'elementwise largest integer <= x',
          TD(flts, f='floor'),
          TD(M, f='floor'),
          ),
'rint' :
    Ufunc(1, 1, None,
          'round x elementwise to the nearest integer, round halfway cases away from zero',
          TD(inexact, f='rint'),
          TD(M, f='rint'),
          ),
'arctan2' :
    Ufunc(2, 1, None,
          'a safe and correct arctan(x1/x2)',
          TD(flts, f='atan2'),
          TD(M, f='arctan2'),
          ),
'remainder' :
    Ufunc(2, 1, None,
          'computes x1-n*x2 where n is floor(x1 / x2)',
          TD(intflt),
          TD(O, f='PyNumber_Remainder'),
          ),
'hypot' :
    Ufunc(2, 1, None,
          'sqrt(x1**2 + x2**2) elementwise',
          TD(flts, f='hypot'),
          TD(M, f='hypot'),
          ),
'isnan' :
    Ufunc(1, 1, None,
          'returns True where x is Not-A-Number',
          TD(inexact, out='?'),
          ),
'isinf' :
    Ufunc(1, 1, None,
          'returns True where x is +inf or -inf',
          TD(inexact, out='?'),
          ),
'isfinite' :
    Ufunc(1, 1, None,
          'returns True where x is finite',
          TD(inexact, out='?'),
          ),
'signbit' :
    Ufunc(1, 1, None,
          'returns True where signbit of x is set (x<0).',
          TD(flts, out='?'),
          ),
'modf' :
    Ufunc(1, 2, None,
          'breaks x into fractional (y1) and integral (y2) parts.\\n\\n    Each output has the same sign as the input.',
          TD(flts),
          ),
}

def indent(st,spaces):
    indention = ' '*spaces
    indented = indention + st.replace('\n','\n'+indention)
    # trim off any trailing spaces
    indented = re.sub(r' +$',r'',indented)
    return indented

chartotype1 = {'f': 'f_f',
               'd': 'd_d',
               'g': 'g_g',
               'F': 'F_F',
               'D': 'D_D',
               'G': 'G_G',
               'O': 'O_O',
               'M': 'O_O_method'}

chartotype2 = {'f': 'ff_f',
               'd': 'dd_d',
               'g': 'gg_g',
               'F': 'FF_F',
               'D': 'DD_D',
               'G': 'GG_G',
               'O': 'OO_O',
               'M': 'OO_O_method'}
#for each name
# 1) create functions, data, and signature
# 2) fill in functions and data in InitOperators
# 3) add function.

# String-handling utilities to avoid locale-dependence.

import string
UPPER_TABLE = string.maketrans(string.ascii_lowercase, string.ascii_uppercase)

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


def make_arrays(funcdict):
    # functions array contains an entry for every type implemented
    #   NULL should be placed where PyUfunc_ style function will be filled in later
    #
    code1list = []
    code2list = []
    names = funcdict.keys()
    names.sort()
    for name in names:
        uf = funcdict[name]
        funclist = []
        datalist = []
        siglist = []
        k = 0
        sub = 0

        if uf.nin > 1:
            assert uf.nin == 2
            thedict = chartotype2  # two inputs and one output
        else:
            thedict = chartotype1  # one input and one output

        for t in uf.type_descriptions:
            if t.func_data is not None:
                funclist.append('NULL')
                astr = '%s_functions[%d] = PyUFunc_%s;' % \
                       (name, k, thedict[t.type])
                code2list.append(astr)
                if t.type == 'O':
                    astr = '%s_data[%d] = (void *) %s;' % \
                           (name, k, t.func_data)
                    code2list.append(astr)
                    datalist.append('(void *)NULL')
                elif t.type == 'M':
                    datalist.append('(void *)"%s"' % t.func_data)
                else:
                    astr = '%s_data[%d] = (void *) %s;' % \
                           (name, k, t.func_data)
                    code2list.append(astr)
                    datalist.append('(void *)NULL')
                    #datalist.append('(void *)%s' % t.func_data)
                sub += 1
            else:
                datalist.append('(void *)NULL');
                tname = english_upper(chartoname[t.type])
                funclist.append('%s_%s' % (tname, name))

            for x in t.in_ + t.out:
                siglist.append('PyArray_%s' % (english_upper(chartoname[x]),))

            k += 1

        funcnames = ', '.join(funclist)
        signames = ', '.join(siglist)
        datanames = ', '.join(datalist)
        code1list.append("static PyUFuncGenericFunction %s_functions[] = { %s };" \
                         % (name, funcnames))
        code1list.append("static void * %s_data[] = { %s };" \
                         % (name, datanames))
        code1list.append("static char %s_signatures[] = { %s };" \
                         % (name, signames))
    return "\n".join(code1list),"\n".join(code2list)

def make_ufuncs(funcdict):
    code3list = []
    names = funcdict.keys()
    names.sort()
    for name in names:
        uf = funcdict[name]
        mlist = []
        mlist.append(\
r"""f = PyUFunc_FromFuncAndData(%s_functions, %s_data, %s_signatures, %d,
                                %d, %d, %s, "%s",
                                "%s", 0);""" % (name, name, name,
                                                len(uf.type_descriptions),
                                                uf.nin, uf.nout,
                                                uf.identity,
                                                name, uf.docstring))
        mlist.append(r"""PyDict_SetItemString(dictionary, "%s", f);""" % name)
        mlist.append(r"""Py_DECREF(f);""")
        code3list.append('\n'.join(mlist))
    return '\n'.join(code3list)


def make_code(funcdict,filename):
    code1, code2 = make_arrays(funcdict)
    code3 = make_ufuncs(funcdict)
    code2 = indent(code2,4)
    code3 = indent(code3,4)
    code = r"""

/** Warning this file is autogenerated!!!

    Please make changes to the code generator program (%s)
**/

%s

static void
InitOperators(PyObject *dictionary) {
    PyObject *f;

%s
%s
}
""" % (filename, code1, code2, code3)
    return code;


if __name__ == "__main__":
    filename = __file__
    fid = open('__umath_generated.c','w')
    code = make_code(defdict, filename)
    fid.write(code)
    fid.close()
