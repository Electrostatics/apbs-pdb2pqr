__docformat__ = "restructuredtext en"
__all__ = ['logspace', 'linspace',
           'select', 'piecewise', 'trim_zeros',
           'copy', 'iterable',
           'diff', 'gradient', 'angle', 'unwrap', 'sort_complex', 'disp',
           'unique', 'extract', 'place', 'nansum', 'nanmax', 'nanargmax',
           'nanargmin', 'nanmin', 'vectorize', 'asarray_chkfinite', 'average',
           'histogram', 'histogramdd', 'bincount', 'digitize', 'cov',
           'corrcoef', 'msort', 'median', 'sinc', 'hamming', 'hanning',
           'bartlett', 'blackman', 'kaiser', 'trapz', 'i0', 'add_newdoc',
           'add_docstring', 'meshgrid', 'delete', 'insert', 'append',
           'interp'
           ]
import warnings

import types
import numpy.core.numeric as _nx
from numpy.core.numeric import ones, zeros, arange, concatenate, array, \
     asarray, asanyarray, empty, empty_like, ndarray, around
from numpy.core.numeric import ScalarType, dot, where, newaxis, intp, \
     integer, isscalar
from numpy.core.umath import pi, multiply, add, arctan2,  \
     frompyfunc, isnan, cos, less_equal, sqrt, sin, mod, exp, log10
from numpy.core.fromnumeric import ravel, nonzero, choose, sort, mean
from numpy.core.numerictypes import typecodes, number
from numpy.lib.shape_base import atleast_1d, atleast_2d
from numpy.lib.twodim_base import diag
from _compiled_base import _insert, add_docstring
from _compiled_base import digitize, bincount, interp as compiled_interp
from arraysetops import setdiff1d
import numpy as np

#end Fernando's utilities

def linspace(start, stop, num=50, endpoint=True, retstep=False):
    """Return evenly spaced numbers.

    Return num evenly spaced samples from start to stop.  If
    endpoint is True, the last sample is stop. If retstep is
    True then return (seq, step_value), where step_value used.

    Parameters
    ----------
    start : {float}
        The value the sequence starts at.
    stop : {float}
        The value the sequence stops at. If ``endpoint`` is false, then
        this is not included in the sequence. Otherwise it is
        guaranteed to be the last value.
    num : {integer}
        Number of samples to generate. Default is 50.
    endpoint : {boolean}
        If true, ``stop`` is the last sample. Otherwise, it is not
        included. Default is true.
    retstep : {boolean}
        If true, return ``(samples, step)``, where ``step`` is the
        spacing used in generating the samples.

    Returns
    -------
    samples : {array}
        ``num`` equally spaced samples from the range [start, stop]
        or [start, stop).
    step : {float} (Only if ``retstep`` is true)
        Size of spacing between samples.

    See Also
    --------
    arange : Similiar to linspace, however, when used with
             a float endpoint, that endpoint may or may not be included.
    logspace

    """
    num = int(num)
    if num <= 0:
        return array([], float)
    if endpoint:
        if num == 1:
            return array([float(start)])
        step = (stop-start)/float((num-1))
        y = _nx.arange(0, num) * step + start
        y[-1] = stop
    else:
        step = (stop-start)/float(num)
        y = _nx.arange(0, num) * step + start
    if retstep:
        return y, step
    else:
        return y

def logspace(start,stop,num=50,endpoint=True,base=10.0):
    """Evenly spaced numbers on a logarithmic scale.

    Computes int(num) evenly spaced exponents from base**start to
    base**stop. If endpoint=True, then last number is base**stop
    """
    y = linspace(start,stop,num=num,endpoint=endpoint)
    return _nx.power(base,y)

def iterable(y):
    try: iter(y)
    except: return 0
    return 1

def histogram(a, bins=10, range=None, normed=False, weights=None, new=False):
    """Compute the histogram from a set of data.

    Parameters
    ----------
    a : array
        The data to histogram.

    bins : int or sequence
        If an int, then the number of equal-width bins in the given
        range.  If new=True, bins can also be the bin edges, allowing
        for non-constant bin widths.

    range : (float, float)
        The lower and upper range of the bins. If not provided, range
        is simply (a.min(), a.max()). Using new=False, lower than
        range are ignored, and values higher than range are tallied in
        the rightmost bin. Using new=True, both lower and upper
        outliers are ignored.

    normed : bool
        If False, the result array will contain the number of samples
        in each bin.  If True, the result array is the value of the
        probability *density* function at the bin normalized such that
        the *integral* over the range is 1. Note that the sum of all
        of the histogram values will not usually be 1; it is not a
        probability *mass* function.

    weights : array
        An array of weights, the same shape as a. If normed is False,
        the histogram is computed by summing the weights of the values
        falling into each bin. If normed is True, the weights are
        normalized, so that the integral of the density over the range
        is 1. This option is only available with new=True.

    new : bool
        Compatibility argument to transition from the old version
        (v1.1) to the new version (v1.2).

    Returns
    -------
    hist : array
        The values of the histogram. See `normed` and `weights` for a
        description of the possible semantics.

    bin_edges : float array
        With new=False, return the left bin edges (length(hist)).
        With new=True, return the bin edges (length(hist)+1).

    See Also
    --------
    histogramdd

    """
    # Old behavior
    if new is False:
        warnings.warn("""
        The semantics of histogram will be modified in
        release 1.2 to improve outlier handling. The new behavior can be
        obtained using new=True. Note that the new version accepts/returns
        the bin edges instead of the left bin edges.
        Please read the docstring for more information.""", FutureWarning)
        a = asarray(a).ravel()

        if (range is not None):
            mn, mx = range
            if (mn > mx):
                raise AttributeError, \
                    'max must be larger than min in range parameter.'

        if not iterable(bins):
            if range is None:
                range = (a.min(), a.max())
            else:
                warnings.warn("""
                Outliers handling will change in version 1.2.
                Please read the docstring for details.""", FutureWarning)
            mn, mx = [mi+0.0 for mi in range]
            if mn == mx:
                mn -= 0.5
                mx += 0.5
            bins = linspace(mn, mx, bins, endpoint=False)
        else:
            if normed:
                raise ValueError, 'Use new=True to pass bin edges explicitly.'
            warnings.warn("""
            The semantic for bins will change in version 1.2.
            The bins will become the bin edges, instead of the left bin edges.
            """, FutureWarning)
            bins = asarray(bins)
            if (np.diff(bins) < 0).any():
                raise AttributeError, 'bins must increase monotonically.'


        if weights is not None:
            raise ValueError, 'weights are only available with new=True.'

        # best block size probably depends on processor cache size
        block = 65536
        n = sort(a[:block]).searchsorted(bins)
        for i in xrange(block, a.size, block):
            n += sort(a[i:i+block]).searchsorted(bins)
        n = concatenate([n, [len(a)]])
        n = n[1:]-n[:-1]

        if normed:
            db = bins[1] - bins[0]
            return 1.0/(a.size*db) * n, bins
        else:
            return n, bins



    # New behavior
    elif new is True:
        a = asarray(a)
        if weights is not None:
            weights = asarray(weights)
            if np.any(weights.shape != a.shape):
                raise ValueError, 'weights should have the same shape as a.'
            weights = weights.ravel()
        a =  a.ravel()

        if (range is not None):
            mn, mx = range
            if (mn > mx):
                raise AttributeError, \
                    'max must be larger than min in range parameter.'

        if not iterable(bins):
            if range is None:
                range = (a.min(), a.max())
            mn, mx = [mi+0.0 for mi in range]
            if mn == mx:
                mn -= 0.5
                mx += 0.5
            bins = linspace(mn, mx, bins+1, endpoint=True)
        else:
            bins = asarray(bins)
            if (np.diff(bins) < 0).any():
                raise AttributeError, 'bins must increase monotonically.'

        # Histogram is an integer or a float array depending on the weights.
        if weights is None:
            ntype = int
        else:
            ntype = weights.dtype
        n = np.zeros(bins.shape, ntype)

        block = 65536
        if weights is None:
            for i in arange(0, len(a), block):
                sa = sort(a[i:i+block])
                n += np.r_[sa.searchsorted(bins[:-1], 'left'), \
                    sa.searchsorted(bins[-1], 'right')]
        else:
            zero = array(0, dtype=ntype)
            for i in arange(0, len(a), block):
                tmp_a = a[i:i+block]
                tmp_w = weights[i:i+block]
                sorting_index = np.argsort(tmp_a)
                sa = tmp_a[sorting_index]
                sw = tmp_w[sorting_index]
                cw = np.concatenate(([zero,], sw.cumsum()))
                bin_index = np.r_[sa.searchsorted(bins[:-1], 'left'), \
                    sa.searchsorted(bins[-1], 'right')]
                n += cw[bin_index]

        n = np.diff(n)

        if normed is False:
            return n, bins
        elif normed is True:
            db = array(np.diff(bins), float)
            return n/(n*db).sum(), bins


def histogramdd(sample, bins=10, range=None, normed=False, weights=None):
    """histogramdd(sample, bins=10, range=None, normed=False, weights=None)

    Return the N-dimensional histogram of the sample.

    Parameters
    ----------
    sample : sequence or array
        A sequence containing N arrays or an NxM array. Input data.

    bins : sequence or scalar
        A sequence of edge arrays, a sequence of bin counts, or a scalar
        which is the bin count for all dimensions. Default is 10.

    range : sequence
        A sequence of lower and upper bin edges. Default is [min, max].

    normed : boolean
        If False, return the number of samples in each bin, if True,
        returns the density.

    weights : array
        Array of weights.  The weights are normed only if normed is True.
        Should the sum of the weights not equal N, the total bin count will
        not be equal to the number of samples.

    Returns
    -------
    hist : array
        Histogram array.

    edges : list
        List of arrays defining the lower bin edges.

    See Also
    --------
    histogram

    Examples
    --------
    >>> x = random.randn(100,3)
    >>> hist3d, edges = histogramdd(x, bins = (5, 6, 7))

    """

    try:
        # Sample is an ND-array.
        N, D = sample.shape
    except (AttributeError, ValueError):
        # Sample is a sequence of 1D arrays.
        sample = atleast_2d(sample).T
        N, D = sample.shape

    nbin = empty(D, int)
    edges = D*[None]
    dedges = D*[None]
    if weights is not None:
        weights = asarray(weights)

    try:
        M = len(bins)
        if M != D:
            raise AttributeError, 'The dimension of bins must be equal ' \
                                  'to the dimension of the sample x.'
    except TypeError:
        bins = D*[bins]

    # Select range for each dimension
    # Used only if number of bins is given.
    if range is None:
        smin = atleast_1d(array(sample.min(0), float))
        smax = atleast_1d(array(sample.max(0), float))
    else:
        smin = zeros(D)
        smax = zeros(D)
        for i in arange(D):
            smin[i], smax[i] = range[i]

    # Make sure the bins have a finite width.
    for i in arange(len(smin)):
        if smin[i] == smax[i]:
            smin[i] = smin[i] - .5
            smax[i] = smax[i] + .5

    # Create edge arrays
    for i in arange(D):
        if isscalar(bins[i]):
            nbin[i] = bins[i] + 2 # +2 for outlier bins
            edges[i] = linspace(smin[i], smax[i], nbin[i]-1)
        else:
            edges[i] = asarray(bins[i], float)
            nbin[i] = len(edges[i])+1  # +1 for outlier bins
        dedges[i] = diff(edges[i])

    nbin =  asarray(nbin)

    # Compute the bin number each sample falls into.
    Ncount = {}
    for i in arange(D):
        Ncount[i] = digitize(sample[:,i], edges[i])

    # Using digitize, values that fall on an edge are put in the right bin.
    # For the rightmost bin, we want values equal to the right
    # edge to be counted in the last bin, and not as an outlier.
    outliers = zeros(N, int)
    for i in arange(D):
        # Rounding precision
        decimal = int(-log10(dedges[i].min())) +6
        # Find which points are on the rightmost edge.
        on_edge = where(around(sample[:,i], decimal) == around(edges[i][-1],
                                                               decimal))[0]
        # Shift these points one bin to the left.
        Ncount[i][on_edge] -= 1

    # Flattened histogram matrix (1D)
    hist = zeros(nbin.prod(), float)

    # Compute the sample indices in the flattened histogram matrix.
    ni = nbin.argsort()
    shape = []
    xy = zeros(N, int)
    for i in arange(0, D-1):
        xy += Ncount[ni[i]] * nbin[ni[i+1:]].prod()
    xy += Ncount[ni[-1]]

    # Compute the number of repetitions in xy and assign it to the
    # flattened histmat.
    if len(xy) == 0:
        return zeros(nbin-2, int), edges

    flatcount = bincount(xy, weights)
    a = arange(len(flatcount))
    hist[a] = flatcount

    # Shape into a proper matrix
    hist = hist.reshape(sort(nbin))
    for i in arange(nbin.size):
        j = ni.argsort()[i]
        hist = hist.swapaxes(i,j)
        ni[i],ni[j] = ni[j],ni[i]

    # Remove outliers (indices 0 and -1 for each dimension).
    core = D*[slice(1,-1)]
    hist = hist[core]

    # Normalize if normed is True
    if normed:
        s = hist.sum()
        for i in arange(D):
            shape = ones(D, int)
            shape[i] = nbin[i]-2
            hist = hist / dedges[i].reshape(shape)
        hist /= s

    if (hist.shape != nbin-2).any():
        raise 'Internal Shape Error'
    return hist, edges


def average(a, axis=None, weights=None, returned=False):
    """Return the weighted average of array a over the given axis.


    Parameters
    ----------
    a : array_like
        Data to be averaged.
    axis : {None, integer}, optional
        Axis along which to average a. If None, averaging is done over the
        entire array irrespective of its shape.
    weights : {None, array_like}, optional
        The importance each datum has in the computation of the
        average. The weights array can either be 1D, in which case  its length
        must be the size of a along the given axis, or of the same shape as a.
        If weights=None, all data are assumed to have weight equal to one.
    returned :{False, boolean}, optional
        If True, the tuple (average, sum_of_weights) is returned,
        otherwise only the average is returmed. Note that if weights=None, then
        the sum of the weights is also the number of elements averaged over.

    Returns
    -------
    average, [sum_of_weights] : {array_type, double}
        Return the average along the specified axis. When returned is True,
        return a tuple with the average as the first element and the sum
        of the weights as the second element. The return type is Float if a is
        of integer type, otherwise it is of the same type as a.
        sum_of_weights is has the same type as the average.


    Examples
    --------
      >>> average(range(1,11), weights=range(10,0,-1))
      4.0

    Raises
    ------
    ZeroDivisionError
        When all weights along axis are zero. See numpy.ma.average for a
        version robust to this type of error.
    TypeError
        When the length of 1D weights is not the same as the shape of a
        along axis.

    """
    if not isinstance(a, np.matrix) :
        a = np.asarray(a)

    if weights is None :
        avg = a.mean(axis)
        scl = avg.dtype.type(a.size/avg.size)
    else :
        a = a + 0.0
        wgt = np.array(weights, dtype=a.dtype, copy=0)

        # Sanity checks
        if a.shape != wgt.shape :
            if axis is None :
                raise TypeError, "Axis must be specified when shapes of a and weights differ."
            if wgt.ndim != 1 :
                raise TypeError, "1D weights expected when shapes of a and weights differ."
            if wgt.shape[0] != a.shape[axis] :
                raise ValueError, "Length of weights not compatible with specified axis."

            # setup wgt to broadcast along axis
            wgt = np.array(wgt, copy=0, ndmin=a.ndim).swapaxes(-1,axis)

        scl = wgt.sum(axis=axis)
        if (scl == 0.0).any():
            raise ZeroDivisionError, "Weights sum to zero, can't be normalized"

        avg = np.multiply(a,wgt).sum(axis)/scl

    if returned:
        scl = np.multiply(avg,0) + scl
        return avg, scl
    else:
        return avg

def asarray_chkfinite(a):
    """Like asarray, but check that no NaNs or Infs are present.
    """
    a = asarray(a)
    if (a.dtype.char in typecodes['AllFloat']) \
           and (_nx.isnan(a).any() or _nx.isinf(a).any()):
        raise ValueError, "array must not contain infs or NaNs"
    return a

def piecewise(x, condlist, funclist, *args, **kw):
    """Return a piecewise-defined function.

    x is the domain

    condlist is a list of boolean arrays or a single boolean array
      The length of the condition list must be n2 or n2-1 where n2
      is the length of the function list.  If len(condlist)==n2-1, then
      an 'otherwise' condition is formed by |'ing all the conditions
      and inverting.

    funclist is a list of functions to call of length (n2).
      Each function should return an array output for an array input
      Each function can take (the same set) of extra arguments and
      keyword arguments which are passed in after the function list.
      A constant may be used in funclist for a function that returns a
      constant (e.g. val  and lambda x: val are equivalent in a funclist).

    The output is the same shape and type as x and is found by
      calling the functions on the appropriate portions of x.

    Note: This is similar to choose or select, except
          the the functions are only evaluated on elements of x
          that satisfy the corresponding condition.

    The result is
           |--
           |  f1(x)  for condition1
     y = --|  f2(x)  for condition2
           |   ...
           |  fn(x)  for conditionn
           |--

    """
    x = asanyarray(x)
    n2 = len(funclist)
    if not isinstance(condlist, type([])):
        condlist = [condlist]
    n = len(condlist)
    if n == n2-1:  # compute the "otherwise" condition.
        totlist = condlist[0]
        for k in range(1, n):
            totlist |= condlist[k]
        condlist.append(~totlist)
        n += 1
    if (n != n2):
        raise ValueError, "function list and condition list must be the same"
    y = empty(x.shape, x.dtype)
    for k in range(n):
        item = funclist[k]
        if not callable(item):
            y[condlist[k]] = item
        else:
            y[condlist[k]] = item(x[condlist[k]], *args, **kw)
    return y

def select(condlist, choicelist, default=0):
    """Return an array composed of different elements in choicelist,
    depending on the list of conditions.

    :Parameters:
        condlist : list of N boolean arrays of length M
            The conditions C_0 through C_(N-1) which determine
            from which vector the output elements are taken.
        choicelist : list of N arrays of length M
            Th vectors V_0 through V_(N-1), from which the output
            elements are chosen.

    :Returns:
        output : 1-dimensional array of length M
            The output at position m is the m-th element of the first
            vector V_n for which C_n[m] is non-zero.  Note that the
            output depends on the order of conditions, since the
            first satisfied condition is used.

            Equivalent to:

                output = []
                for m in range(M):
                    output += [V[m] for V,C in zip(values,cond) if C[m]]
                              or [default]

    """
    n = len(condlist)
    n2 = len(choicelist)
    if n2 != n:
        raise ValueError, "list of cases must be same length as list of conditions"
    choicelist = [default] + choicelist
    S = 0
    pfac = 1
    for k in range(1, n+1):
        S += k * pfac * asarray(condlist[k-1])
        if k < n:
            pfac *= (1-asarray(condlist[k-1]))
    # handle special case of a 1-element condition but
    #  a multi-element choice
    if type(S) in ScalarType or max(asarray(S).shape)==1:
        pfac = asarray(1)
        for k in range(n2+1):
            pfac = pfac + asarray(choicelist[k])
        if type(S) in ScalarType:
            S = S*ones(asarray(pfac).shape, type(S))
        else:
            S = S*ones(asarray(pfac).shape, S.dtype)
    return choose(S, tuple(choicelist))

def _asarray1d(arr, copy=False):
    """Ensure 1D array for one array.
    """
    if copy:
        return asarray(arr).flatten()
    else:
        return asarray(arr).ravel()

def copy(a):
    """Return an array copy of the given object.
    """
    return array(a, copy=True)

# Basic operations

def gradient(f, *varargs):
    """Calculate the gradient of an N-dimensional scalar function.

    Uses central differences on the interior and first differences on boundaries
    to give the same shape.

    Inputs:

      f -- An N-dimensional array giving samples of a scalar function

      varargs -- 0, 1, or N scalars giving the sample distances in each direction

    Outputs:

      N arrays of the same shape as f giving the derivative of f with respect
      to each dimension.

    """
    N = len(f.shape)  # number of dimensions
    n = len(varargs)
    if n == 0:
        dx = [1.0]*N
    elif n == 1:
        dx = [varargs[0]]*N
    elif n == N:
        dx = list(varargs)
    else:
        raise SyntaxError, "invalid number of arguments"

    # use central differences on interior and first differences on endpoints

    outvals = []

    # create slice objects --- initially all are [:, :, ..., :]
    slice1 = [slice(None)]*N
    slice2 = [slice(None)]*N
    slice3 = [slice(None)]*N

    otype = f.dtype.char
    if otype not in ['f', 'd', 'F', 'D']:
        otype = 'd'

    for axis in range(N):
        # select out appropriate parts for this dimension
        out = zeros(f.shape, f.dtype.char)
        slice1[axis] = slice(1, -1)
        slice2[axis] = slice(2, None)
        slice3[axis] = slice(None, -2)
        # 1D equivalent -- out[1:-1] = (f[2:] - f[:-2])/2.0
        out[slice1] = (f[slice2] - f[slice3])/2.0
        slice1[axis] = 0
        slice2[axis] = 1
        slice3[axis] = 0
        # 1D equivalent -- out[0] = (f[1] - f[0])
        out[slice1] = (f[slice2] - f[slice3])
        slice1[axis] = -1
        slice2[axis] = -1
        slice3[axis] = -2
        # 1D equivalent -- out[-1] = (f[-1] - f[-2])
        out[slice1] = (f[slice2] - f[slice3])

        # divide by step size
        outvals.append(out / dx[axis])

        # reset the slice object in this dimension to ":"
        slice1[axis] = slice(None)
        slice2[axis] = slice(None)
        slice3[axis] = slice(None)

    if N == 1:
        return outvals[0]
    else:
        return outvals


def diff(a, n=1, axis=-1):
    """Calculate the nth order discrete difference along given axis.
    """
    if n == 0:
        return a
    if n < 0:
        raise ValueError, 'order must be non-negative but got ' + repr(n)
    a = asanyarray(a)
    nd = len(a.shape)
    slice1 = [slice(None)]*nd
    slice2 = [slice(None)]*nd
    slice1[axis] = slice(1, None)
    slice2[axis] = slice(None, -1)
    slice1 = tuple(slice1)
    slice2 = tuple(slice2)
    if n > 1:
        return diff(a[slice1]-a[slice2], n-1, axis=axis)
    else:
        return a[slice1]-a[slice2]

try:
    add_docstring(digitize,
r"""digitize(x,bins)

Return the index of the bin to which each value of x belongs.

Each index i returned is such that bins[i-1] <= x < bins[i] if
bins is monotonically increasing, or bins [i-1] > x >= bins[i] if
bins is monotonically decreasing.

Beyond the bounds of the bins 0 or len(bins) is returned as appropriate.

""")
except RuntimeError:
    pass

try:
    add_docstring(bincount,
r"""bincount(x,weights=None)

Return the number of occurrences of each value in x.

x must be a list of non-negative integers.  The output, b[i],
represents the number of times that i is found in x.  If weights
is specified, every occurrence of i at a position p contributes
weights[p] instead of 1.

See also: histogram, digitize, unique.

""")
except RuntimeError:
    pass

try:
    add_docstring(add_docstring,
r"""docstring(obj, docstring)

Add a docstring to a built-in obj if possible.
If the obj already has a docstring raise a RuntimeError
If this routine does not know how to add a docstring to the object
raise a TypeError

""")
except RuntimeError:
    pass


def interp(x, xp, fp, left=None, right=None):
    """Return the value of a piecewise-linear function at each value in x.

    The piecewise-linear function, f, is defined by the known data-points
    fp=f(xp). The xp points must be sorted in increasing order but this is
    not checked.

    For values of x < xp[0] return the value given by left.  If left is None,
    then return fp[0].
    For values of x > xp[-1] return the value given by right. If right is
    None, then return fp[-1].
    """
    if isinstance(x, (float, int, number)):
        return compiled_interp([x], xp, fp, left, right).item()
    else:
        return compiled_interp(x, xp, fp, left, right)


def angle(z, deg=0):
    """
    Return the angle of the complex argument z.

    Examples
    --------
    >>> numpy.angle(1+1j)          # in radians
    0.78539816339744828
    >>> numpy.angle(1+1j,deg=True) # in degrees
    45.0

    """
    if deg:
        fact = 180/pi
    else:
        fact = 1.0
    z = asarray(z)
    if (issubclass(z.dtype.type, _nx.complexfloating)):
        zimag = z.imag
        zreal = z.real
    else:
        zimag = 0
        zreal = z
    return arctan2(zimag, zreal) * fact

def unwrap(p, discont=pi, axis=-1):
    """Unwrap radian phase p by changing absolute jumps greater than
       'discont' to their 2*pi complement along the given axis.
    """
    p = asarray(p)
    nd = len(p.shape)
    dd = diff(p, axis=axis)
    slice1 = [slice(None, None)]*nd     # full slices
    slice1[axis] = slice(1, None)
    ddmod = mod(dd+pi, 2*pi)-pi
    _nx.putmask(ddmod, (ddmod==-pi) & (dd > 0), pi)
    ph_correct = ddmod - dd;
    _nx.putmask(ph_correct, abs(dd)<discont, 0)
    up = array(p, copy=True, dtype='d')
    up[slice1] = p[slice1] + ph_correct.cumsum(axis)
    return up

def sort_complex(a):
    """ Sort 'a' as a complex array using the real part first and then
    the imaginary part if the real part is equal (the default sort order
    for complex arrays).  This function is a wrapper ensuring a complex
    return type.

    """
    b = array(a,copy=True)
    b.sort()
    if not issubclass(b.dtype.type, _nx.complexfloating):
        if b.dtype.char in 'bhBH':
            return b.astype('F')
        elif b.dtype.char == 'g':
            return b.astype('G')
        else:
            return b.astype('D')
    else:
        return b

def trim_zeros(filt, trim='fb'):
    """ Trim the leading and trailing zeros from a 1D array.

    Examples
    --------
    >>> import numpy
    >>> a = array((0, 0, 0, 1, 2, 3, 2, 1, 0))
    >>> numpy.trim_zeros(a)
    array([1, 2, 3, 2, 1])

    """
    first = 0
    trim = trim.upper()
    if 'F' in trim:
        for i in filt:
            if i != 0.: break
            else: first = first + 1
    last = len(filt)
    if 'B' in trim:
        for i in filt[::-1]:
            if i != 0.: break
            else: last = last - 1
    return filt[first:last]

import sys
if sys.hexversion < 0x2040000:
    from sets import Set as set

def unique(x):
    """
    Return sorted unique items from an array or sequence.

    Examples
    --------
    >>> numpy.unique([5,2,4,0,4,4,2,2,1])
    array([0, 1, 2, 4, 5])

    """
    try:
        tmp = x.flatten()
        if tmp.size == 0:
            return tmp
        tmp.sort()
        idx = concatenate(([True],tmp[1:]!=tmp[:-1]))
        return tmp[idx]
    except AttributeError:
        items = list(set(x))
        items.sort()
        return asarray(items)

def extract(condition, arr):
    """Return the elements of ravel(arr) where ravel(condition) is True
    (in 1D).

    Equivalent to compress(ravel(condition), ravel(arr)).
    """
    return _nx.take(ravel(arr), nonzero(ravel(condition))[0])

def place(arr, mask, vals):
    """Similar to putmask arr[mask] = vals but the 1D array vals has the
    same number of elements as the non-zero values of mask. Inverse of
    extract.

    """
    return _insert(arr, mask, vals)

def nansum(a, axis=None):
    """Sum the array over the given axis, treating NaNs as 0.
    """
    y = array(a,subok=True)
    if not issubclass(y.dtype.type, _nx.integer):
        y[isnan(a)] = 0
    return y.sum(axis)

def nanmin(a, axis=None):
    """Find the minimium over the given axis, ignoring NaNs.
    """
    y = array(a,subok=True)
    if not issubclass(y.dtype.type, _nx.integer):
        y[isnan(a)] = _nx.inf
    return y.min(axis)

def nanargmin(a, axis=None):
    """Find the indices of the minimium over the given axis ignoring NaNs.
    """
    y = array(a, subok=True)
    if not issubclass(y.dtype.type, _nx.integer):
        y[isnan(a)] = _nx.inf
    return y.argmin(axis)

def nanmax(a, axis=None):
    """Find the maximum over the given axis ignoring NaNs.
    """
    y = array(a, subok=True)
    if not issubclass(y.dtype.type, _nx.integer):
        y[isnan(a)] = -_nx.inf
    return y.max(axis)

def nanargmax(a, axis=None):
    """Find the maximum over the given axis ignoring NaNs.
    """
    y = array(a,subok=True)
    if not issubclass(y.dtype.type, _nx.integer):
        y[isnan(a)] = -_nx.inf
    return y.argmax(axis)

def disp(mesg, device=None, linefeed=True):
    """Display a message to the given device (default is sys.stdout)
    with or without a linefeed.
    """
    if device is None:
        import sys
        device = sys.stdout
    if linefeed:
        device.write('%s\n' % mesg)
    else:
        device.write('%s' % mesg)
    device.flush()
    return

# return number of input arguments and
#  number of default arguments
import re
def _get_nargs(obj):
    if not callable(obj):
        raise TypeError, "Object is not callable."
    if hasattr(obj,'func_code'):
        fcode = obj.func_code
        nargs = fcode.co_argcount
        if obj.func_defaults is not None:
            ndefaults = len(obj.func_defaults)
        else:
            ndefaults = 0
        if isinstance(obj, types.MethodType):
            nargs -= 1
        return nargs, ndefaults
    terr = re.compile(r'.*? takes exactly (?P<exargs>\d+) argument(s|) \((?P<gargs>\d+) given\)')
    try:
        obj()
        return 0, 0
    except TypeError, msg:
        m = terr.match(str(msg))
        if m:
            nargs = int(m.group('exargs'))
            ndefaults = int(m.group('gargs'))
            if isinstance(obj, types.MethodType):
                nargs -= 1
            return nargs, ndefaults
    raise ValueError, 'failed to determine the number of arguments for %s' % (obj)


class vectorize(object):
    """
    vectorize(somefunction, otypes=None, doc=None)

    Generalized function class.

    Define a vectorized function which takes nested sequence
    of objects or numpy arrays as inputs and returns a
    numpy array as output, evaluating the function over successive
    tuples of the input arrays like the python map function except it uses
    the broadcasting rules of numpy.

    Data-type of output of vectorized is determined by calling the function
    with the first element of the input.  This can be avoided by specifying
    the otypes argument as either a string of typecode characters or a list
    of data-types specifiers.  There should be one data-type specifier for
    each output.

    Parameters
    ----------
    f : callable
      A Python function or method.

    Examples
    --------
    >>> def myfunc(a, b):
    ...    if a > b:
    ...        return a-b
    ...    else:
    ...        return a+b

    >>> vfunc = vectorize(myfunc)

    >>> vfunc([1, 2, 3, 4], 2)
    array([3, 4, 1, 2])

    """
    def __init__(self, pyfunc, otypes='', doc=None):
        self.thefunc = pyfunc
        self.ufunc = None
        nin, ndefault = _get_nargs(pyfunc)
        if nin == 0 and ndefault == 0:
            self.nin = None
            self.nin_wo_defaults = None
        else:
            self.nin = nin
            self.nin_wo_defaults = nin - ndefault
        self.nout = None
        if doc is None:
            self.__doc__ = pyfunc.__doc__
        else:
            self.__doc__ = doc
        if isinstance(otypes, types.StringType):
            self.otypes = otypes
            for char in self.otypes:
                if char not in typecodes['All']:
                    raise ValueError, "invalid otype specified"
        elif iterable(otypes):
            self.otypes = ''.join([_nx.dtype(x).char for x in otypes])
        else:
            raise ValueError, "output types must be a string of typecode characters or a list of data-types"
        self.lastcallargs = 0

    def __call__(self, *args):
        # get number of outputs and output types by calling
        #  the function on the first entries of args
        nargs = len(args)
        if self.nin:
            if (nargs > self.nin) or (nargs < self.nin_wo_defaults):
                raise ValueError, "mismatch between python function inputs"\
                      " and received arguments"

        # we need a new ufunc if this is being called with more arguments.
        if (self.lastcallargs != nargs):
            self.lastcallargs = nargs
            self.ufunc = None
            self.nout = None

        if self.nout is None or self.otypes == '':
            newargs = []
            for arg in args:
                newargs.append(asarray(arg).flat[0])
            theout = self.thefunc(*newargs)
            if isinstance(theout, types.TupleType):
                self.nout = len(theout)
            else:
                self.nout = 1
                theout = (theout,)
            if self.otypes == '':
                otypes = []
                for k in range(self.nout):
                    otypes.append(asarray(theout[k]).dtype.char)
                self.otypes = ''.join(otypes)

        # Create ufunc if not already created
        if (self.ufunc is None):
            self.ufunc = frompyfunc(self.thefunc, nargs, self.nout)

        # Convert to object arrays first
        newargs = [array(arg,copy=False,subok=True,dtype=object) for arg in args]
        if self.nout == 1:
            _res = array(self.ufunc(*newargs),copy=False,
                         subok=True,dtype=self.otypes[0])
        else:
            _res = tuple([array(x,copy=False,subok=True,dtype=c) \
                          for x, c in zip(self.ufunc(*newargs), self.otypes)])
        return _res

def cov(m, y=None, rowvar=1, bias=0):
    """Estimate the covariance matrix.

    If m is a vector, return the variance.  For matrices return the
    covariance matrix.

    If y is given it is treated as an additional (set of)
    variable(s).

    Normalization is by (N-1) where N is the number of observations
    (unbiased estimate).  If bias is 1 then normalization is by N.

    If rowvar is non-zero (default), then each row is a variable with
    observations in the columns, otherwise each column
    is a variable and the observations are in the rows.
    """

    X = array(m, ndmin=2, dtype=float)
    if X.shape[0] == 1:
        rowvar = 1
    if rowvar:
        axis = 0
        tup = (slice(None),newaxis)
    else:
        axis = 1
        tup = (newaxis, slice(None))


    if y is not None:
        y = array(y, copy=False, ndmin=2, dtype=float)
        X = concatenate((X,y),axis)

    X -= X.mean(axis=1-axis)[tup]
    if rowvar:
        N = X.shape[1]
    else:
        N = X.shape[0]

    if bias:
        fact = N*1.0
    else:
        fact = N-1.0

    if not rowvar:
        return (dot(X.T, X.conj()) / fact).squeeze()
    else:
        return (dot(X, X.T.conj()) / fact).squeeze()

def corrcoef(x, y=None, rowvar=1, bias=0):
    """The correlation coefficients
    """
    c = cov(x, y, rowvar, bias)
    try:
        d = diag(c)
    except ValueError: # scalar covariance
        return 1
    return c/sqrt(multiply.outer(d,d))

def blackman(M):
    """blackman(M) returns the M-point Blackman window.
    """
    if M < 1:
        return array([])
    if M == 1:
        return ones(1, float)
    n = arange(0,M)
    return 0.42-0.5*cos(2.0*pi*n/(M-1)) + 0.08*cos(4.0*pi*n/(M-1))

def bartlett(M):
    """
    Return the Bartlett window.

    The Bartlett window is very similar to a triangular window, except
    that the end points are at zero.  It is often used in signal
    processing for tapering a signal, without generating too much
    ripple in the frequency domain.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an
        empty array is returned.

    Returns
    -------
    out : array
        The triangular window, normalized to one (the value one
        appears only if the number of samples is odd), with the first
        and last samples equal to zero.

    See Also
    --------
    blackman, hamming, hanning, kaiser

    Notes
    -----
    The Bartlett window is defined as

    .. math:: w(n) = \\frac{2}{M-1} \left(
              \\frac{M-1}{2} - \\left|n - \\frac{M-1}{2}\\right|
              \\right)

    Most references to the Bartlett window come from the signal
    processing literature, where it is used as one of many windowing
    functions for smoothing values.  Note that convolution with this
    window produces linear interpolation.  It is also known as an
    apodization (which means"removing the foot", i.e. smoothing
    discontinuities at the beginning and end of the sampled signal) or
    tapering function.

    References
    ----------
    .. [1] M.S. Bartlett, "Periodogram Analysis and Continuous Spectra",
           Biometrika 37, 1-16, 1950.
    .. [2] A.V. Oppenheim and R.W. Schafer, "Discrete-Time Signal
           Processing", Prentice-Hall, 1999, pp. 468-471.
    .. [3] Wikipedia, "Window function",
           http://en.wikipedia.org/wiki/Window_function
    .. [4] W.H. Press,  B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
           "Numerical Recipes", Cambridge University Press, 1986, page 429.

    Examples
    --------
    >>> from numpy import bartlett
    >>> bartlett(12)
    array([ 0.        ,  0.18181818,  0.36363636,  0.54545455,  0.72727273,
            0.90909091,  0.90909091,  0.72727273,  0.54545455,  0.36363636,
            0.18181818,  0.        ])

    Plot the window and its frequency response:

    >>> from numpy import clip, log10, array, bartlett
    >>> from scipy.fftpack import fft
    >>> from matplotlib import pyplot as plt

    >>> window = bartlett(51)
    >>> plt.plot(window)
    >>> plt.title("Bartlett window")
    >>> plt.ylabel("Amplitude")
    >>> plt.xlabel("Sample")
    >>> plt.show()

    >>> A = fft(window, 2048) / 25.5
    >>> mag = abs(fftshift(A))
    >>> freq = linspace(-0.5,0.5,len(A))
    >>> response = 20*log10(mag)
    >>> response = clip(response,-100,100)
    >>> plt.plot(freq, response)
    >>> plt.title("Frequency response of Bartlett window")
    >>> plt.ylabel("Magnitude [dB]")
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    >>> plt.axis('tight'); plt.show()

    """
    if M < 1:
        return array([])
    if M == 1:
        return ones(1, float)
    n = arange(0,M)
    return where(less_equal(n,(M-1)/2.0),2.0*n/(M-1),2.0-2.0*n/(M-1))

def hanning(M):
    """hanning(M) returns the M-point Hanning window.
    """
    if M < 1:
        return array([])
    if M == 1:
        return ones(1, float)
    n = arange(0,M)
    return 0.5-0.5*cos(2.0*pi*n/(M-1))

def hamming(M):
    """hamming(M) returns the M-point Hamming window.
    """
    if M < 1:
        return array([])
    if M == 1:
        return ones(1,float)
    n = arange(0,M)
    return 0.54-0.46*cos(2.0*pi*n/(M-1))

## Code from cephes for i0

_i0A = [
-4.41534164647933937950E-18,
 3.33079451882223809783E-17,
-2.43127984654795469359E-16,
 1.71539128555513303061E-15,
-1.16853328779934516808E-14,
 7.67618549860493561688E-14,
-4.85644678311192946090E-13,
 2.95505266312963983461E-12,
-1.72682629144155570723E-11,
 9.67580903537323691224E-11,
-5.18979560163526290666E-10,
 2.65982372468238665035E-9,
-1.30002500998624804212E-8,
 6.04699502254191894932E-8,
-2.67079385394061173391E-7,
 1.11738753912010371815E-6,
-4.41673835845875056359E-6,
 1.64484480707288970893E-5,
-5.75419501008210370398E-5,
 1.88502885095841655729E-4,
-5.76375574538582365885E-4,
 1.63947561694133579842E-3,
-4.32430999505057594430E-3,
 1.05464603945949983183E-2,
-2.37374148058994688156E-2,
 4.93052842396707084878E-2,
-9.49010970480476444210E-2,
 1.71620901522208775349E-1,
-3.04682672343198398683E-1,
 6.76795274409476084995E-1]

_i0B = [
-7.23318048787475395456E-18,
-4.83050448594418207126E-18,
 4.46562142029675999901E-17,
 3.46122286769746109310E-17,
-2.82762398051658348494E-16,
-3.42548561967721913462E-16,
 1.77256013305652638360E-15,
 3.81168066935262242075E-15,
-9.55484669882830764870E-15,
-4.15056934728722208663E-14,
 1.54008621752140982691E-14,
 3.85277838274214270114E-13,
 7.18012445138366623367E-13,
-1.79417853150680611778E-12,
-1.32158118404477131188E-11,
-3.14991652796324136454E-11,
 1.18891471078464383424E-11,
 4.94060238822496958910E-10,
 3.39623202570838634515E-9,
 2.26666899049817806459E-8,
 2.04891858946906374183E-7,
 2.89137052083475648297E-6,
 6.88975834691682398426E-5,
 3.36911647825569408990E-3,
 8.04490411014108831608E-1]

def _chbevl(x, vals):
    b0 = vals[0]
    b1 = 0.0

    for i in xrange(1,len(vals)):
        b2 = b1
        b1 = b0
        b0 = x*b1 - b2 + vals[i]

    return 0.5*(b0 - b2)

def _i0_1(x):
    return exp(x) * _chbevl(x/2.0-2, _i0A)

def _i0_2(x):
    return exp(x) * _chbevl(32.0/x - 2.0, _i0B) / sqrt(x)

def i0(x):
    x = atleast_1d(x).copy()
    y = empty_like(x)
    ind = (x<0)
    x[ind] = -x[ind]
    ind = (x<=8.0)
    y[ind] = _i0_1(x[ind])
    ind2 = ~ind
    y[ind2] = _i0_2(x[ind2])
    return y.squeeze()

## End of cephes code for i0

def kaiser(M,beta):
    """kaiser(M, beta) returns a Kaiser window of length M with shape parameter
    beta.
    """
    from numpy.dual import i0
    n = arange(0,M)
    alpha = (M-1)/2.0
    return i0(beta * sqrt(1-((n-alpha)/alpha)**2.0))/i0(beta)

def sinc(x):
    """sinc(x) returns sin(pi*x)/(pi*x) at all points of array x.
    """
    y = pi* where(x == 0, 1.0e-20, x)
    return sin(y)/y

def msort(a):
    b = array(a,subok=True,copy=True)
    b.sort(0)
    return b

def median(a, axis=0, out=None, overwrite_input=False):
    """Compute the median along the specified axis.

    Returns the median of the array elements.  The median is taken
    over the first axis of the array by default, otherwise over
    the specified axis.

    Parameters
    ----------
    a : array-like
        Input array or object that can be converted to an array
    axis : {int, None}, optional
        Axis along which the medians are computed. The default is to
        compute the median along the first dimension.  axis=None
        returns the median of the flattened array

    out : ndarray, optional
        Alternative output array in which to place the result. It must
        have the same shape and buffer length as the expected output
        but the type will be cast if necessary.

    overwrite_input : {False, True}, optional
       If True, then allow use of memory of input array (a) for
       calculations. The input array will be modified by the call to
       median. This will save memory when you do not need to preserve
       the contents of the input array. Treat the input as undefined,
       but it will probably be fully or partially sorted. Default is
       False. Note that, if overwrite_input is true, and the input
       is not already an ndarray, an error will be raised.

    Returns
    -------
    median : ndarray.
        A new array holding the result is returned unless out is
        specified, in which case a reference to out is returned.
        Return datatype is float64 for ints and floats smaller than
        float64, or the input datatype otherwise.

    See Also
    -------
    mean

    Notes
    -----
    Given a vector V length N, the median of V is the middle value of
    a sorted copy of V (Vs) - i.e. Vs[(N-1)/2], when N is odd. It is
    the mean of the two middle values of Vs, when N is even.

    Examples
    --------
    >>> import numpy as np
    >>> from numpy import median
    >>> a = np.array([[10, 7, 4], [3, 2, 1]])
    >>> a
    array([[10,  7,  4],
           [ 3,  2,  1]])
    >>> median(a)
    array([ 6.5,  4.5,  2.5])
    >>> median(a, axis=None)
    3.5
    >>> median(a, axis=1)
    array([ 7.,  2.])
    >>> m = median(a)
    >>> out = np.zeros_like(m)
    >>> median(a, out=m)
    array([ 6.5,  4.5,  2.5])
    >>> m
    array([ 6.5,  4.5,  2.5])
    >>> b = a.copy()
    >>> median(b, axis=1, overwrite_input=True)
    array([ 7.,  2.])
    >>> assert not np.all(a==b)
    >>> b = a.copy()
    >>> median(b, axis=None, overwrite_input=True)
    3.5
    >>> assert not np.all(a==b)
    """
    if overwrite_input:
        if axis is None:
            sorted = a.ravel()
            sorted.sort()
        else:
            a.sort(axis=axis)
            sorted = a
    else:
        sorted = sort(a, axis=axis)
    if axis is None:
        axis = 0
    indexer = [slice(None)] * sorted.ndim
    index = int(sorted.shape[axis]/2)
    if sorted.shape[axis] % 2 == 1:
        # index with slice to allow mean (below) to work
        indexer[axis] = slice(index, index+1)
    else:
        indexer[axis] = slice(index-1, index+1)
    # Use mean in odd and even case to coerce data type
    # and check, use out array.
    return mean(sorted[indexer], axis=axis, out=out)

def trapz(y, x=None, dx=1.0, axis=-1):
    """Integrate y(x) using samples along the given axis and the composite
    trapezoidal rule.  If x is None, spacing given by dx is assumed.
    """
    y = asarray(y)
    if x is None:
        d = dx
    else:
        d = diff(x,axis=axis)
    nd = len(y.shape)
    slice1 = [slice(None)]*nd
    slice2 = [slice(None)]*nd
    slice1[axis] = slice(1,None)
    slice2[axis] = slice(None,-1)
    return add.reduce(d * (y[slice1]+y[slice2])/2.0,axis)

#always succeed
def add_newdoc(place, obj, doc):
    """Adds documentation to obj which is in module place.

    If doc is a string add it to obj as a docstring

    If doc is a tuple, then the first element is interpreted as
       an attribute of obj and the second as the docstring
          (method, docstring)

    If doc is a list, then each element of the list should be a
       sequence of length two --> [(method1, docstring1),
       (method2, docstring2), ...]

    This routine never raises an error.
       """
    try:
        new = {}
        exec 'from %s import %s' % (place, obj) in new
        if isinstance(doc, str):
            add_docstring(new[obj], doc.strip())
        elif isinstance(doc, tuple):
            add_docstring(getattr(new[obj], doc[0]), doc[1].strip())
        elif isinstance(doc, list):
            for val in doc:
                add_docstring(getattr(new[obj], val[0]), val[1].strip())
    except:
        pass


# From matplotlib
def meshgrid(x,y):
    """
    For vectors x, y with lengths Nx=len(x) and Ny=len(y), return X, Y
    where X and Y are (Ny, Nx) shaped arrays with the elements of x
    and y repeated to fill the matrix

    EG,

      [X, Y] = meshgrid([1,2,3], [4,5,6,7])

       X =
         1   2   3
         1   2   3
         1   2   3
         1   2   3


       Y =
         4   4   4
         5   5   5
         6   6   6
         7   7   7
  """
    x = asarray(x)
    y = asarray(y)
    numRows, numCols = len(y), len(x)  # yes, reversed
    x = x.reshape(1,numCols)
    X = x.repeat(numRows, axis=0)

    y = y.reshape(numRows,1)
    Y = y.repeat(numCols, axis=1)
    return X, Y

def delete(arr, obj, axis=None):
    """Return a new array with sub-arrays along an axis deleted.

    Return a new array with the sub-arrays (i.e. rows or columns)
    deleted along the given axis as specified by obj

    obj may be a slice_object (s_[3:5:2]) or an integer
    or an array of integers indicated which sub-arrays to
    remove.

    If axis is None, then ravel the array first.

    Examples
    --------
    >>> arr = [[3,4,5],
    ...       [1,2,3],
    ...       [6,7,8]]

    >>> delete(arr, 1, 1)
    array([[3, 5],
           [1, 3],
           [6, 8]])
    >>> delete(arr, 1, 0)
    array([[3, 4, 5],
           [6, 7, 8]])
    """
    wrap = None
    if type(arr) is not ndarray:
        try:
            wrap = arr.__array_wrap__
        except AttributeError:
            pass


    arr = asarray(arr)
    ndim = arr.ndim
    if axis is None:
        if ndim != 1:
            arr = arr.ravel()
        ndim = arr.ndim;
        axis = ndim-1;
    if ndim == 0:
        if wrap:
            return wrap(arr)
        else:
            return arr.copy()
    slobj = [slice(None)]*ndim
    N = arr.shape[axis]
    newshape = list(arr.shape)
    if isinstance(obj, (int, long, integer)):
        if (obj < 0): obj += N
        if (obj < 0 or obj >=N):
            raise ValueError, "invalid entry"
        newshape[axis]-=1;
        new = empty(newshape, arr.dtype, arr.flags.fnc)
        slobj[axis] = slice(None, obj)
        new[slobj] = arr[slobj]
        slobj[axis] = slice(obj,None)
        slobj2 = [slice(None)]*ndim
        slobj2[axis] = slice(obj+1,None)
        new[slobj] = arr[slobj2]
    elif isinstance(obj, slice):
        start, stop, step = obj.indices(N)
        numtodel = len(xrange(start, stop, step))
        if numtodel <= 0:
            if wrap:
                return wrap(new)
            else:
                return arr.copy()
        newshape[axis] -= numtodel
        new = empty(newshape, arr.dtype, arr.flags.fnc)
        # copy initial chunk
        if start == 0:
            pass
        else:
            slobj[axis] = slice(None, start)
            new[slobj] = arr[slobj]
        # copy end chunck
        if stop == N:
            pass
        else:
            slobj[axis] = slice(stop-numtodel,None)
            slobj2 = [slice(None)]*ndim
            slobj2[axis] = slice(stop, None)
            new[slobj] = arr[slobj2]
        # copy middle pieces
        if step == 1:
            pass
        else:  # use array indexing.
            obj = arange(start, stop, step, dtype=intp)
            all = arange(start, stop, dtype=intp)
            obj = setdiff1d(all, obj)
            slobj[axis] = slice(start, stop-numtodel)
            slobj2 = [slice(None)]*ndim
            slobj2[axis] = obj
            new[slobj] = arr[slobj2]
    else: # default behavior
        obj = array(obj, dtype=intp, copy=0, ndmin=1)
        all = arange(N, dtype=intp)
        obj = setdiff1d(all, obj)
        slobj[axis] = obj
        new = arr[slobj]
    if wrap:
        return wrap(new)
    else:
        return new

def insert(arr, obj, values, axis=None):
    """Return a new array with values inserted along the given axis
    before the given indices

    If axis is None, then ravel the array first.

    The obj argument can be an integer, a slice, or a sequence of
    integers.

    Examples
    --------
    >>> a = array([[1,2,3],
    ...            [4,5,6],
    ...            [7,8,9]])

    >>> insert(a, [1,2], [[4],[5]], axis=0)
    array([[1, 2, 3],
           [4, 4, 4],
           [4, 5, 6],
           [5, 5, 5],
           [7, 8, 9]])
    """
    wrap = None
    if type(arr) is not ndarray:
        try:
            wrap = arr.__array_wrap__
        except AttributeError:
            pass

    arr = asarray(arr)
    ndim = arr.ndim
    if axis is None:
        if ndim != 1:
            arr = arr.ravel()
        ndim = arr.ndim
        axis = ndim-1
    if (ndim == 0):
        arr = arr.copy()
        arr[...] = values
        if wrap:
            return wrap(arr)
        else:
            return arr
    slobj = [slice(None)]*ndim
    N = arr.shape[axis]
    newshape = list(arr.shape)
    if isinstance(obj, (int, long, integer)):
        if (obj < 0): obj += N
        if obj < 0 or obj > N:
            raise ValueError, "index (%d) out of range (0<=index<=%d) "\
                  "in dimension %d" % (obj, N, axis)
        newshape[axis] += 1;
        new = empty(newshape, arr.dtype, arr.flags.fnc)
        slobj[axis] = slice(None, obj)
        new[slobj] = arr[slobj]
        slobj[axis] = obj
        new[slobj] = values
        slobj[axis] = slice(obj+1,None)
        slobj2 = [slice(None)]*ndim
        slobj2[axis] = slice(obj,None)
        new[slobj] = arr[slobj2]
        if wrap:
            return wrap(new)
        return new

    elif isinstance(obj, slice):
        # turn it into a range object
        obj = arange(*obj.indices(N),**{'dtype':intp})

    # get two sets of indices
    #  one is the indices which will hold the new stuff
    #  two is the indices where arr will be copied over

    obj = asarray(obj, dtype=intp)
    numnew = len(obj)
    index1 = obj + arange(numnew)
    index2 = setdiff1d(arange(numnew+N),index1)
    newshape[axis] += numnew
    new = empty(newshape, arr.dtype, arr.flags.fnc)
    slobj2 = [slice(None)]*ndim
    slobj[axis] = index1
    slobj2[axis] = index2
    new[slobj] = values
    new[slobj2] = arr

    if wrap:
        return wrap(new)
    return new

def append(arr, values, axis=None):
    """Append to the end of an array along axis (ravel first if None)
    """
    arr = asanyarray(arr)
    if axis is None:
        if arr.ndim != 1:
            arr = arr.ravel()
        values = ravel(values)
        axis = arr.ndim-1
    return concatenate((arr, values), axis=axis)
