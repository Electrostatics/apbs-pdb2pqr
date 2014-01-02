# pylint: disable-msg=W0611, W0612, W0511
"""Tests suite for MaskedArray.
Adapted from the original test_ma by Pierre Gerard-Marchant

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id: test_extras.py 3473 2007-10-29 15:18:13Z jarrod.millman $
"""
__author__ = "Pierre GF Gerard-Marchant ($Author: jarrod.millman $)"
__version__ = '1.0'
__revision__ = "$Revision: 3473 $"
__date__     = '$Date: 2007-10-29 17:18:13 +0200 (Mon, 29 Oct 2007) $'

import numpy as N
from numpy.testing import NumpyTest, NumpyTestCase
from numpy.testing.utils import build_err_msg

import numpy.ma.testutils
from numpy.ma.testutils import *

import numpy.ma.core
from numpy.ma.core import *
import numpy.ma.extras
from numpy.ma.extras import *

class TestAverage(NumpyTestCase):
    "Several tests of average. Why so many ? Good point..."
    def check_testAverage1(self):
        "Test of average."
        ott = array([0.,1.,2.,3.], mask=[1,0,0,0])
        assert_equal(2.0, average(ott,axis=0))
        assert_equal(2.0, average(ott, weights=[1., 1., 2., 1.]))
        result, wts = average(ott, weights=[1.,1.,2.,1.], returned=1)
        assert_equal(2.0, result)
        assert(wts == 4.0)
        ott[:] = masked
        assert_equal(average(ott,axis=0).mask, [True])
        ott = array([0.,1.,2.,3.], mask=[1,0,0,0])
        ott = ott.reshape(2,2)
        ott[:,1] = masked
        assert_equal(average(ott,axis=0), [2.0, 0.0])
        assert_equal(average(ott,axis=1).mask[0], [True])
        assert_equal([2.,0.], average(ott, axis=0))
        result, wts = average(ott, axis=0, returned=1)
        assert_equal(wts, [1., 0.])

    def check_testAverage2(self):
        "More tests of average."
        w1 = [0,1,1,1,1,0]
        w2 = [[0,1,1,1,1,0],[1,0,0,0,0,1]]
        x = arange(6, dtype=float_)
        assert_equal(average(x, axis=0), 2.5)
        assert_equal(average(x, axis=0, weights=w1), 2.5)
        y = array([arange(6, dtype=float_), 2.0*arange(6)])
        assert_equal(average(y, None), N.add.reduce(N.arange(6))*3./12.)
        assert_equal(average(y, axis=0), N.arange(6) * 3./2.)
        assert_equal(average(y, axis=1), [average(x,axis=0), average(x,axis=0) * 2.0])
        assert_equal(average(y, None, weights=w2), 20./6.)
        assert_equal(average(y, axis=0, weights=w2), [0.,1.,2.,3.,4.,10.])
        assert_equal(average(y, axis=1), [average(x,axis=0), average(x,axis=0) * 2.0])
        m1 = zeros(6)
        m2 = [0,0,1,1,0,0]
        m3 = [[0,0,1,1,0,0],[0,1,1,1,1,0]]
        m4 = ones(6)
        m5 = [0, 1, 1, 1, 1, 1]
        assert_equal(average(masked_array(x, m1),axis=0), 2.5)
        assert_equal(average(masked_array(x, m2),axis=0), 2.5)
        assert_equal(average(masked_array(x, m4),axis=0).mask, [True])
        assert_equal(average(masked_array(x, m5),axis=0), 0.0)
        assert_equal(count(average(masked_array(x, m4),axis=0)), 0)
        z = masked_array(y, m3)
        assert_equal(average(z, None), 20./6.)
        assert_equal(average(z, axis=0), [0.,1.,99.,99.,4.0, 7.5])
        assert_equal(average(z, axis=1), [2.5, 5.0])
        assert_equal(average(z,axis=0, weights=w2), [0.,1., 99., 99., 4.0, 10.0])

    def check_testAverage3(self):
        "Yet more tests of average!"
        a = arange(6)
        b = arange(6) * 3
        r1, w1 = average([[a,b],[b,a]], axis=1, returned=1)
        assert_equal(shape(r1) , shape(w1))
        assert_equal(r1.shape , w1.shape)
        r2, w2 = average(ones((2,2,3)), axis=0, weights=[3,1], returned=1)
        assert_equal(shape(w2) , shape(r2))
        r2, w2 = average(ones((2,2,3)), returned=1)
        assert_equal(shape(w2) , shape(r2))
        r2, w2 = average(ones((2,2,3)), weights=ones((2,2,3)), returned=1)
        assert_equal(shape(w2), shape(r2))
        a2d = array([[1,2],[0,4]], float)
        a2dm = masked_array(a2d, [[0,0],[1,0]])
        a2da = average(a2d, axis=0)
        assert_equal(a2da, [0.5, 3.0])
        a2dma = average(a2dm, axis=0)
        assert_equal(a2dma, [1.0, 3.0])
        a2dma = average(a2dm, axis=None)
        assert_equal(a2dma, 7./3.)
        a2dma = average(a2dm, axis=1)
        assert_equal(a2dma, [1.5, 4.0])

class TestConcatenator(NumpyTestCase):
    "Tests for mr_, the equivalent of r_ for masked arrays."
    def check_1d(self):
        "Tests mr_ on 1D arrays."
        assert_array_equal(mr_[1,2,3,4,5,6],array([1,2,3,4,5,6]))
        b = ones(5)
        m = [1,0,0,0,0]
        d = masked_array(b,mask=m)
        c = mr_[d,0,0,d]
        assert(isinstance(c,MaskedArray) or isinstance(c,core.MaskedArray))
        assert_array_equal(c,[1,1,1,1,1,0,0,1,1,1,1,1])
        assert_array_equal(c.mask, mr_[m,0,0,m])

    def check_2d(self):
        "Tests mr_ on 2D arrays."
        a_1 = rand(5,5)
        a_2 = rand(5,5)
        m_1 = N.round_(rand(5,5),0)
        m_2 = N.round_(rand(5,5),0)
        b_1 = masked_array(a_1,mask=m_1)
        b_2 = masked_array(a_2,mask=m_2)
        d = mr_['1',b_1,b_2]  # append columns
        assert(d.shape == (5,10))
        assert_array_equal(d[:,:5],b_1)
        assert_array_equal(d[:,5:],b_2)
        assert_array_equal(d.mask, N.r_['1',m_1,m_2])
        d = mr_[b_1,b_2]
        assert(d.shape == (10,5))
        assert_array_equal(d[:5,:],b_1)
        assert_array_equal(d[5:,:],b_2)
        assert_array_equal(d.mask, N.r_[m_1,m_2])

class TestNotMasked(NumpyTestCase):
    "Tests notmasked_edges and notmasked_contiguous."
    def check_edges(self):
        "Tests unmasked_edges"
        a = masked_array(N.arange(24).reshape(3,8),
                         mask=[[0,0,0,0,1,1,1,0],
                               [1,1,1,1,1,1,1,1],
                               [0,0,0,0,0,0,1,0],])
        #
        assert_equal(notmasked_edges(a, None), [0,23])
        #
        tmp = notmasked_edges(a, 0)
        assert_equal(tmp[0], (array([0,0,0,0,2,2,0]), array([0,1,2,3,4,5,7])))
        assert_equal(tmp[1], (array([2,2,2,2,2,2,2]), array([0,1,2,3,4,5,7])))
        #
        tmp = notmasked_edges(a, 1)
        assert_equal(tmp[0], (array([0,2,]), array([0,0])))
        assert_equal(tmp[1], (array([0,2,]), array([7,7])))

    def check_contiguous(self):
        "Tests notmasked_contiguous"
        a = masked_array(N.arange(24).reshape(3,8),
                         mask=[[0,0,0,0,1,1,1,1],
                               [1,1,1,1,1,1,1,1],
                               [0,0,0,0,0,0,1,0],])
        tmp = notmasked_contiguous(a, None)
        assert_equal(tmp[-1], slice(23,23,None))
        assert_equal(tmp[-2], slice(16,21,None))
        assert_equal(tmp[-3], slice(0,3,None))
        #
        tmp = notmasked_contiguous(a, 0)
        assert(len(tmp[-1]) == 1)
        assert(tmp[-2] is None)
        assert_equal(tmp[-3],tmp[-1])
        assert(len(tmp[0]) == 2)
        #
        tmp = notmasked_contiguous(a, 1)
        assert_equal(tmp[0][-1], slice(0,3,None))
        assert(tmp[1] is None)
        assert_equal(tmp[2][-1], slice(7,7,None))
        assert_equal(tmp[2][-2], slice(0,5,None))

class Test2DFunctions(NumpyTestCase):
    "Tests 2D functions"
    def check_compress2d(self):
        "Tests compress2d"
        x = array(N.arange(9).reshape(3,3), mask=[[1,0,0],[0,0,0],[0,0,0]])
        assert_equal(compress_rowcols(x), [[4,5],[7,8]] )
        assert_equal(compress_rowcols(x,0), [[3,4,5],[6,7,8]] )
        assert_equal(compress_rowcols(x,1), [[1,2],[4,5],[7,8]] )
        x = array(x._data, mask=[[0,0,0],[0,1,0],[0,0,0]])
        assert_equal(compress_rowcols(x), [[0,2],[6,8]] )
        assert_equal(compress_rowcols(x,0), [[0,1,2],[6,7,8]] )
        assert_equal(compress_rowcols(x,1), [[0,2],[3,5],[6,8]] )
        x = array(x._data, mask=[[1,0,0],[0,1,0],[0,0,0]])
        assert_equal(compress_rowcols(x), [[8]] )
        assert_equal(compress_rowcols(x,0), [[6,7,8]] )
        assert_equal(compress_rowcols(x,1,), [[2],[5],[8]] )
        x = array(x._data, mask=[[1,0,0],[0,1,0],[0,0,1]])
        assert_equal(compress_rowcols(x).size, 0 )
        assert_equal(compress_rowcols(x,0).size, 0 )
        assert_equal(compress_rowcols(x,1).size, 0 )
    #
    def check_mask_rowcols(self):
        "Tests mask_rowcols."
        x = array(N.arange(9).reshape(3,3), mask=[[1,0,0],[0,0,0],[0,0,0]])
        assert_equal(mask_rowcols(x).mask, [[1,1,1],[1,0,0],[1,0,0]] )
        assert_equal(mask_rowcols(x,0).mask, [[1,1,1],[0,0,0],[0,0,0]] )
        assert_equal(mask_rowcols(x,1).mask, [[1,0,0],[1,0,0],[1,0,0]] )
        x = array(x._data, mask=[[0,0,0],[0,1,0],[0,0,0]])
        assert_equal(mask_rowcols(x).mask, [[0,1,0],[1,1,1],[0,1,0]] )
        assert_equal(mask_rowcols(x,0).mask, [[0,0,0],[1,1,1],[0,0,0]] )
        assert_equal(mask_rowcols(x,1).mask, [[0,1,0],[0,1,0],[0,1,0]] )
        x = array(x._data, mask=[[1,0,0],[0,1,0],[0,0,0]])
        assert_equal(mask_rowcols(x).mask, [[1,1,1],[1,1,1],[1,1,0]] )
        assert_equal(mask_rowcols(x,0).mask, [[1,1,1],[1,1,1],[0,0,0]] )
        assert_equal(mask_rowcols(x,1,).mask, [[1,1,0],[1,1,0],[1,1,0]] )
        x = array(x._data, mask=[[1,0,0],[0,1,0],[0,0,1]])
        assert(mask_rowcols(x).all())
        assert(mask_rowcols(x,0).all())
        assert(mask_rowcols(x,1).all())
    #
    def test_dot(self):
        "Tests dot product"
        n = N.arange(1,7)
        #
        m = [1,0,0,0,0,0]
        a = masked_array(n, mask=m).reshape(2,3)
        b = masked_array(n, mask=m).reshape(3,2)
        c = dot(a,b,True)
        assert_equal(c.mask, [[1,1],[1,0]])
        c = dot(b,a,True)
        assert_equal(c.mask, [[1,1,1],[1,0,0],[1,0,0]])
        c = dot(a,b,False)
        assert_equal(c, N.dot(a.filled(0), b.filled(0)))
        c = dot(b,a,False)
        assert_equal(c, N.dot(b.filled(0), a.filled(0)))
        #
        m = [0,0,0,0,0,1]
        a = masked_array(n, mask=m).reshape(2,3)
        b = masked_array(n, mask=m).reshape(3,2)
        c = dot(a,b,True)
        assert_equal(c.mask,[[0,1],[1,1]])
        c = dot(b,a,True)
        assert_equal(c.mask, [[0,0,1],[0,0,1],[1,1,1]])
        c = dot(a,b,False)
        assert_equal(c, N.dot(a.filled(0), b.filled(0)))
        assert_equal(c, dot(a,b))
        c = dot(b,a,False)
        assert_equal(c, N.dot(b.filled(0), a.filled(0)))
        #
        m = [0,0,0,0,0,0]
        a = masked_array(n, mask=m).reshape(2,3)
        b = masked_array(n, mask=m).reshape(3,2)
        c = dot(a,b)
        assert_equal(c.mask,nomask)
        c = dot(b,a)
        assert_equal(c.mask,nomask)
        #
        a = masked_array(n, mask=[1,0,0,0,0,0]).reshape(2,3)
        b = masked_array(n, mask=[0,0,0,0,0,0]).reshape(3,2)
        c = dot(a,b,True)
        assert_equal(c.mask,[[1,1],[0,0]])
        c = dot(a,b,False)
        assert_equal(c, N.dot(a.filled(0),b.filled(0)))
        c = dot(b,a,True)
        assert_equal(c.mask,[[1,0,0],[1,0,0],[1,0,0]])
        c = dot(b,a,False)
        assert_equal(c, N.dot(b.filled(0),a.filled(0)))
        #
        a = masked_array(n, mask=[0,0,0,0,0,1]).reshape(2,3)
        b = masked_array(n, mask=[0,0,0,0,0,0]).reshape(3,2)
        c = dot(a,b,True)
        assert_equal(c.mask,[[0,0],[1,1]])
        c = dot(a,b)
        assert_equal(c, N.dot(a.filled(0),b.filled(0)))
        c = dot(b,a,True)
        assert_equal(c.mask,[[0,0,1],[0,0,1],[0,0,1]])
        c = dot(b,a,False)
        assert_equal(c, N.dot(b.filled(0), a.filled(0)))
        #
        a = masked_array(n, mask=[0,0,0,0,0,1]).reshape(2,3)
        b = masked_array(n, mask=[0,0,1,0,0,0]).reshape(3,2)
        c = dot(a,b,True)
        assert_equal(c.mask,[[1,0],[1,1]])
        c = dot(a,b,False)
        assert_equal(c, N.dot(a.filled(0),b.filled(0)))
        c = dot(b,a,True)
        assert_equal(c.mask,[[0,0,1],[1,1,1],[0,0,1]])
        c = dot(b,a,False)
        assert_equal(c, N.dot(b.filled(0),a.filled(0)))

    def test_ediff1d(self):
        "Tests mediff1d"
        x = masked_array(N.arange(5), mask=[1,0,0,0,1])
        difx_d = (x._data[1:]-x._data[:-1])
        difx_m = (x._mask[1:]-x._mask[:-1])
        dx = ediff1d(x)
        assert_equal(dx._data, difx_d)
        assert_equal(dx._mask, difx_m)
        #
        dx = ediff1d(x, to_begin=masked)
        assert_equal(dx._data, N.r_[0,difx_d])
        assert_equal(dx._mask, N.r_[1,difx_m])
        dx = ediff1d(x, to_begin=[1,2,3])
        assert_equal(dx._data, N.r_[[1,2,3],difx_d])
        assert_equal(dx._mask, N.r_[[0,0,0],difx_m])
        #
        dx = ediff1d(x, to_end=masked)
        assert_equal(dx._data, N.r_[difx_d,0])
        assert_equal(dx._mask, N.r_[difx_m,1])
        dx = ediff1d(x, to_end=[1,2,3])
        assert_equal(dx._data, N.r_[difx_d,[1,2,3]])
        assert_equal(dx._mask, N.r_[difx_m,[0,0,0]])
        #
        dx = ediff1d(x, to_end=masked, to_begin=masked)
        assert_equal(dx._data, N.r_[0,difx_d,0])
        assert_equal(dx._mask, N.r_[1,difx_m,1])
        dx = ediff1d(x, to_end=[1,2,3], to_begin=masked)
        assert_equal(dx._data, N.r_[0,difx_d,[1,2,3]])
        assert_equal(dx._mask, N.r_[1,difx_m,[0,0,0]])
        #
        dx = ediff1d(x._data, to_end=masked, to_begin=masked)
        assert_equal(dx._data, N.r_[0,difx_d,0])
        assert_equal(dx._mask, N.r_[1,0,0,0,0,1])

class TestApplyAlongAxis(NumpyTestCase):
    "Tests 2D functions"
    def check_3d(self):
        a = arange(12.).reshape(2,2,3)
        def myfunc(b):
            return b[1]
        xa = apply_along_axis(myfunc,2,a)
        assert_equal(xa,[[1,4],[7,10]])

class TestMedian(NumpyTestCase):
    def __init__(self, *args, **kwds):
        NumpyTestCase.__init__(self, *args, **kwds)
    #
    def test_2d(self):
        "Tests median w/ 2D"
        (n,p) = (101,30)
        x = masked_array(numpy.linspace(-1.,1.,n),)
        x[:10] = x[-10:] = masked
        z = masked_array(numpy.empty((n,p), dtype=numpy.float_))
        z[:,0] = x[:]
        idx = numpy.arange(len(x))
        for i in range(1,p):
            numpy.random.shuffle(idx)
            z[:,i] = x[idx]
        assert_equal(median(z[:,0]), 0)
        assert_equal(median(z), numpy.zeros((p,)))
    #
    def test_3d(self):
        "Tests median w/ 3D"
        x = numpy.ma.arange(24).reshape(3,4,2)
        x[x%3==0] = masked
        assert_equal(median(x,0), [[12,9],[6,15],[12,9],[18,15]])
        x.shape = (4,3,2)
        assert_equal(median(x,0),[[99,10],[11,99],[13,14]])
        x = numpy.ma.arange(24).reshape(4,3,2)
        x[x%5==0] = masked
        assert_equal(median(x,0), [[12,10],[8,9],[16,17]])


class TestPolynomial(NumpyTestCase):
    #
    def test_polyfit(self):
        "Tests polyfit"
        # On ndarrays
        x = numpy.random.rand(10)
        y = numpy.random.rand(20).reshape(-1,2)
        assert_almost_equal(polyfit(x,y,3),numpy.polyfit(x,y,3))
        # ON 1D maskedarrays
        x = x.view(MaskedArray)
        x[0] = masked
        y = y.view(MaskedArray)
        y[0,0] = y[-1,-1] = masked
        #
        (C,R,K,S,D) = polyfit(x,y[:,0],3,full=True)
        (c,r,k,s,d) = numpy.polyfit(x[1:], y[1:,0].compressed(), 3, full=True)
        for (a,a_) in zip((C,R,K,S,D),(c,r,k,s,d)):
            assert_almost_equal(a, a_)
        #
        (C,R,K,S,D) = polyfit(x,y[:,-1],3,full=True)
        (c,r,k,s,d) = numpy.polyfit(x[1:-1], y[1:-1,-1], 3, full=True)
        for (a,a_) in zip((C,R,K,S,D),(c,r,k,s,d)):
            assert_almost_equal(a, a_)
        #
        (C,R,K,S,D) = polyfit(x,y,3,full=True)
        (c,r,k,s,d) = numpy.polyfit(x[1:-1], y[1:-1,:], 3, full=True)
        for (a,a_) in zip((C,R,K,S,D),(c,r,k,s,d)):
            assert_almost_equal(a, a_)


###############################################################################
#------------------------------------------------------------------------------
if __name__ == "__main__":
    NumpyTest().run()
