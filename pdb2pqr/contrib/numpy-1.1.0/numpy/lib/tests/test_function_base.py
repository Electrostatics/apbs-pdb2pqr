import sys

from numpy.testing import *
set_package_path()
import numpy.lib;reload(numpy.lib)
from numpy.lib import *
from numpy.core import *

del sys.path[0]

class TestAny(NumpyTestCase):
    def check_basic(self):
        y1 = [0,0,1,0]
        y2 = [0,0,0,0]
        y3 = [1,0,1,0]
        assert(any(y1))
        assert(any(y3))
        assert(not any(y2))

    def check_nd(self):
        y1 = [[0,0,0],[0,1,0],[1,1,0]]
        assert(any(y1))
        assert_array_equal(sometrue(y1,axis=0),[1,1,0])
        assert_array_equal(sometrue(y1,axis=1),[0,1,1])

class TestAll(NumpyTestCase):
    def check_basic(self):
        y1 = [0,1,1,0]
        y2 = [0,0,0,0]
        y3 = [1,1,1,1]
        assert(not all(y1))
        assert(all(y3))
        assert(not all(y2))
        assert(all(~array(y2)))

    def check_nd(self):
        y1 = [[0,0,1],[0,1,1],[1,1,1]]
        assert(not all(y1))
        assert_array_equal(alltrue(y1,axis=0),[0,0,1])
        assert_array_equal(alltrue(y1,axis=1),[0,0,1])

class TestAverage(NumpyTestCase):
    def check_basic(self):
        y1 = array([1,2,3])
        assert(average(y1,axis=0) == 2.)
        y2 = array([1.,2.,3.])
        assert(average(y2,axis=0) == 2.)
        y3 = [0.,0.,0.]
        assert(average(y3,axis=0) == 0.)

        y4 = ones((4,4))
        y4[0,1] = 0
        y4[1,0] = 2
        assert_almost_equal(y4.mean(0), average(y4, 0))
        assert_almost_equal(y4.mean(1), average(y4, 1))

        y5 = rand(5,5)
        assert_almost_equal(y5.mean(0), average(y5, 0))
        assert_almost_equal(y5.mean(1), average(y5, 1))

        y6 = matrix(rand(5,5))
        assert_array_equal(y6.mean(0), average(y6,0))

    def check_weights(self):
        y = arange(10)
        w = arange(10)
        assert_almost_equal(average(y, weights=w), (arange(10)**2).sum()*1./arange(10).sum())

        y1 = array([[1,2,3],[4,5,6]])
        w0 = [1,2]
        actual = average(y1,weights=w0,axis=0)
        desired = array([3.,4.,5.])
        assert_almost_equal(actual, desired)


        w1 = [0,0,1]
        desired = array([3., 6.])
        assert_almost_equal(average(y1, weights=w1, axis=1), desired)

        # This should raise an error. Can we test for that ?
        # assert_equal(average(y1, weights=w1), 9./2.)


        # 2D Case
        w2 = [[0,0,1],[0,0,2]]
        desired = array([3., 6.])
        assert_array_equal(average(y1, weights=w2, axis=1), desired)

        assert_equal(average(y1, weights=w2), 5.)


    def check_returned(self):
        y = array([[1,2,3],[4,5,6]])

        # No weights
        avg, scl = average(y, returned=True)
        assert_equal(scl, 6.)

        avg, scl = average(y, 0, returned=True)
        assert_array_equal(scl, array([2.,2.,2.]))

        avg, scl = average(y, 1, returned=True)
        assert_array_equal(scl, array([3.,3.]))

        # With weights
        w0 = [1,2]
        avg, scl = average(y, weights=w0, axis=0, returned=True)
        assert_array_equal(scl, array([3., 3., 3.]))

        w1 = [1,2,3]
        avg, scl = average(y, weights=w1, axis=1, returned=True)
        assert_array_equal(scl, array([6., 6.]))

        w2 = [[0,0,1],[1,2,3]]
        avg, scl = average(y, weights=w2, axis=1, returned=True)
        assert_array_equal(scl, array([1.,6.]))


class TestSelect(NumpyTestCase):
    def _select(self,cond,values,default=0):
        output = []
        for m in range(len(cond)):
            output += [V[m] for V,C in zip(values,cond) if C[m]] or [default]
        return output

    def check_basic(self):
        choices = [array([1,2,3]),
                   array([4,5,6]),
                   array([7,8,9])]
        conditions = [array([0,0,0]),
                      array([0,1,0]),
                      array([0,0,1])]
        assert_array_equal(select(conditions,choices,default=15),
                           self._select(conditions,choices,default=15))

        assert_equal(len(choices),3)
        assert_equal(len(conditions),3)

class TestLogspace(NumpyTestCase):
    def check_basic(self):
        y = logspace(0,6)
        assert(len(y)==50)
        y = logspace(0,6,num=100)
        assert(y[-1] == 10**6)
        y = logspace(0,6,endpoint=0)
        assert(y[-1] < 10**6)
        y = logspace(0,6,num=7)
        assert_array_equal(y,[1,10,100,1e3,1e4,1e5,1e6])

class TestLinspace(NumpyTestCase):
    def check_basic(self):
        y = linspace(0,10)
        assert(len(y)==50)
        y = linspace(2,10,num=100)
        assert(y[-1] == 10)
        y = linspace(2,10,endpoint=0)
        assert(y[-1] < 10)
        y,st = linspace(2,10,retstep=1)
        assert_almost_equal(st,8/49.0)
        assert_array_almost_equal(y,mgrid[2:10:50j],13)

    def check_corner(self):
        y = list(linspace(0,1,1))
        assert y == [0.0], y
        y = list(linspace(0,1,2.5))
        assert y == [0.0, 1.0]

    def check_type(self):
        t1 = linspace(0,1,0).dtype
        t2 = linspace(0,1,1).dtype
        t3 = linspace(0,1,2).dtype
        assert_equal(t1, t2)
        assert_equal(t2, t3)

class TestInsert(NumpyTestCase):
    def check_basic(self):
        a = [1,2,3]
        assert_equal(insert(a,0,1), [1,1,2,3])
        assert_equal(insert(a,3,1), [1,2,3,1])
        assert_equal(insert(a,[1,1,1],[1,2,3]), [1,1,2,3,2,3])

class TestAmax(NumpyTestCase):
    def check_basic(self):
        a = [3,4,5,10,-3,-5,6.0]
        assert_equal(amax(a),10.0)
        b = [[3,6.0, 9.0],
             [4,10.0,5.0],
             [8,3.0,2.0]]
        assert_equal(amax(b,axis=0),[8.0,10.0,9.0])
        assert_equal(amax(b,axis=1),[9.0,10.0,8.0])

class TestAmin(NumpyTestCase):
    def check_basic(self):
        a = [3,4,5,10,-3,-5,6.0]
        assert_equal(amin(a),-5.0)
        b = [[3,6.0, 9.0],
             [4,10.0,5.0],
             [8,3.0,2.0]]
        assert_equal(amin(b,axis=0),[3.0,3.0,2.0])
        assert_equal(amin(b,axis=1),[3.0,4.0,2.0])

class TestPtp(NumpyTestCase):
    def check_basic(self):
        a = [3,4,5,10,-3,-5,6.0]
        assert_equal(ptp(a,axis=0),15.0)
        b = [[3,6.0, 9.0],
             [4,10.0,5.0],
             [8,3.0,2.0]]
        assert_equal(ptp(b,axis=0),[5.0,7.0,7.0])
        assert_equal(ptp(b,axis=-1),[6.0,6.0,6.0])

class TestCumsum(NumpyTestCase):
    def check_basic(self):
        ba = [1,2,10,11,6,5,4]
        ba2 = [[1,2,3,4],[5,6,7,9],[10,3,4,5]]
        for ctype in [int8,uint8,int16,uint16,int32,uint32,
                      float32,float64,complex64,complex128]:
            a = array(ba,ctype)
            a2 = array(ba2,ctype)
            assert_array_equal(cumsum(a,axis=0), array([1,3,13,24,30,35,39],ctype))
            assert_array_equal(cumsum(a2,axis=0), array([[1,2,3,4],[6,8,10,13],
                                                         [16,11,14,18]],ctype))
            assert_array_equal(cumsum(a2,axis=1),
                               array([[1,3,6,10],
                                      [5,11,18,27],
                                      [10,13,17,22]],ctype))

class TestProd(NumpyTestCase):
    def check_basic(self):
        ba = [1,2,10,11,6,5,4]
        ba2 = [[1,2,3,4],[5,6,7,9],[10,3,4,5]]
        for ctype in [int16,uint16,int32,uint32,
                      float32,float64,complex64,complex128]:
            a = array(ba,ctype)
            a2 = array(ba2,ctype)
            if ctype in ['1', 'b']:
                self.failUnlessRaises(ArithmeticError, prod, a)
                self.failUnlessRaises(ArithmeticError, prod, a2, 1)
                self.failUnlessRaises(ArithmeticError, prod, a)
            else:
                assert_equal(prod(a,axis=0),26400)
                assert_array_equal(prod(a2,axis=0),
                                   array([50,36,84,180],ctype))
                assert_array_equal(prod(a2,axis=-1),array([24, 1890, 600],ctype))

class TestCumprod(NumpyTestCase):
    def check_basic(self):
        ba = [1,2,10,11,6,5,4]
        ba2 = [[1,2,3,4],[5,6,7,9],[10,3,4,5]]
        for ctype in [int16,uint16,int32,uint32,
                      float32,float64,complex64,complex128]:
            a = array(ba,ctype)
            a2 = array(ba2,ctype)
            if ctype in ['1', 'b']:
                self.failUnlessRaises(ArithmeticError, cumprod, a)
                self.failUnlessRaises(ArithmeticError, cumprod, a2, 1)
                self.failUnlessRaises(ArithmeticError, cumprod, a)
            else:
                assert_array_equal(cumprod(a,axis=-1),
                                   array([1, 2, 20, 220,
                                          1320, 6600, 26400],ctype))
                assert_array_equal(cumprod(a2,axis=0),
                                   array([[ 1,  2,  3,   4],
                                          [ 5, 12, 21,  36],
                                          [50, 36, 84, 180]],ctype))
                assert_array_equal(cumprod(a2,axis=-1),
                                   array([[ 1,  2,   6,   24],
                                          [ 5, 30, 210, 1890],
                                          [10, 30, 120,  600]],ctype))

class TestDiff(NumpyTestCase):
    def check_basic(self):
        x = [1,4,6,7,12]
        out = array([3,2,1,5])
        out2 = array([-1,-1,4])
        out3 = array([0,5])
        assert_array_equal(diff(x),out)
        assert_array_equal(diff(x,n=2),out2)
        assert_array_equal(diff(x,n=3),out3)

    def check_nd(self):
        x = 20*rand(10,20,30)
        out1 = x[:,:,1:] - x[:,:,:-1]
        out2 = out1[:,:,1:] - out1[:,:,:-1]
        out3 = x[1:,:,:] - x[:-1,:,:]
        out4 = out3[1:,:,:] - out3[:-1,:,:]
        assert_array_equal(diff(x),out1)
        assert_array_equal(diff(x,n=2),out2)
        assert_array_equal(diff(x,axis=0),out3)
        assert_array_equal(diff(x,n=2,axis=0),out4)

class TestAngle(NumpyTestCase):
    def check_basic(self):
        x = [1+3j,sqrt(2)/2.0+1j*sqrt(2)/2,1,1j,-1,-1j,1-3j,-1+3j]
        y = angle(x)
        yo = [arctan(3.0/1.0),arctan(1.0),0,pi/2,pi,-pi/2.0,
              -arctan(3.0/1.0),pi-arctan(3.0/1.0)]
        z = angle(x,deg=1)
        zo = array(yo)*180/pi
        assert_array_almost_equal(y,yo,11)
        assert_array_almost_equal(z,zo,11)

class TestTrimZeros(NumpyTestCase):
    """ only testing for integer splits.
    """
    def check_basic(self):
        a= array([0,0,1,2,3,4,0])
        res = trim_zeros(a)
        assert_array_equal(res,array([1,2,3,4]))
    def check_leading_skip(self):
        a= array([0,0,1,0,2,3,4,0])
        res = trim_zeros(a)
        assert_array_equal(res,array([1,0,2,3,4]))
    def check_trailing_skip(self):
        a= array([0,0,1,0,2,3,0,4,0])
        res = trim_zeros(a)
        assert_array_equal(res,array([1,0,2,3,0,4]))


class TestExtins(NumpyTestCase):
    def check_basic(self):
        a = array([1,3,2,1,2,3,3])
        b = extract(a>1,a)
        assert_array_equal(b,[3,2,2,3,3])
    def check_place(self):
        a = array([1,4,3,2,5,8,7])
        place(a,[0,1,0,1,0,1,0],[2,4,6])
        assert_array_equal(a,[1,2,3,4,5,6,7])
    def check_both(self):
        a = rand(10)
        mask = a > 0.5
        ac = a.copy()
        c = extract(mask, a)
        place(a,mask,0)
        place(a,mask,c)
        assert_array_equal(a,ac)

class TestVectorize(NumpyTestCase):
    def check_simple(self):
        def addsubtract(a,b):
            if a > b:
                return a - b
            else:
                return a + b
        f = vectorize(addsubtract)
        r = f([0,3,6,9],[1,3,5,7])
        assert_array_equal(r,[1,6,1,2])
    def check_scalar(self):
        def addsubtract(a,b):
            if a > b:
                return a - b
            else:
                return a + b
        f = vectorize(addsubtract)
        r = f([0,3,6,9],5)
        assert_array_equal(r,[5,8,1,4])
    def check_large(self):
        x = linspace(-3,2,10000)
        f = vectorize(lambda x: x)
        y = f(x)
        assert_array_equal(y, x)

class TestDigitize(NumpyTestCase):
    def check_forward(self):
        x = arange(-6,5)
        bins = arange(-5,5)
        assert_array_equal(digitize(x,bins),arange(11))

    def check_reverse(self):
        x = arange(5,-6,-1)
        bins = arange(5,-5,-1)
        assert_array_equal(digitize(x,bins),arange(11))

    def check_random(self):
        x = rand(10)
        bin = linspace(x.min(), x.max(), 10)
        assert all(digitize(x,bin) != 0)

class TestUnwrap(NumpyTestCase):
    def check_simple(self):
                #check that unwrap removes jumps greather that 2*pi
        assert_array_equal(unwrap([1,1+2*pi]),[1,1])
        #check that unwrap maintans continuity
        assert(all(diff(unwrap(rand(10)*100))<pi))


class TestFilterwindows(NumpyTestCase):
    def check_hanning(self):
        #check symmetry
        w=hanning(10)
        assert_array_almost_equal(w,flipud(w),7)
        #check known value
        assert_almost_equal(sum(w,axis=0),4.500,4)

    def check_hamming(self):
        #check symmetry
        w=hamming(10)
        assert_array_almost_equal(w,flipud(w),7)
        #check known value
        assert_almost_equal(sum(w,axis=0),4.9400,4)

    def check_bartlett(self):
        #check symmetry
        w=bartlett(10)
        assert_array_almost_equal(w,flipud(w),7)
        #check known value
        assert_almost_equal(sum(w,axis=0),4.4444,4)

    def check_blackman(self):
        #check symmetry
        w=blackman(10)
        assert_array_almost_equal(w,flipud(w),7)
        #check known value
        assert_almost_equal(sum(w,axis=0),3.7800,4)


class TestTrapz(NumpyTestCase):
    def check_simple(self):
        r=trapz(exp(-1.0/2*(arange(-10,10,.1))**2)/sqrt(2*pi),dx=0.1)
        #check integral of normal equals 1
        assert_almost_equal(sum(r,axis=0),1,7)

class TestSinc(NumpyTestCase):
    def check_simple(self):
        assert(sinc(0)==1)
        w=sinc(linspace(-1,1,100))
        #check symmetry
        assert_array_almost_equal(w,flipud(w),7)

class TestHistogram(NumpyTestCase):
    import warnings
    warnings.simplefilter('ignore', FutureWarning)
    def check_simple(self):
        n=100
        v=rand(n)
        (a,b)=histogram(v)
        #check if the sum of the bins equals the number of samples
        assert_equal(sum(a,axis=0), n)
        #check that the bin counts are evenly spaced when the data is from a
        # linear function
        (a,b)=histogram(linspace(0,10,100))
        assert_array_equal(a, 10)

    def check_simple_new(self):
        n=100
        v=rand(n)
        (a,b)=histogram(v, new=True)
        #check if the sum of the bins equals the number of samples
        assert_equal(sum(a,axis=0), n)
        #check that the bin counts are evenly spaced when the data is from a
        # linear function
        (a,b)=histogram(linspace(0,10,100), new=True)
        assert_array_equal(a, 10)

    def check_normed_new(self):
        # Check that the integral of the density equals 1.
        n = 100
        v = rand(n)
        a,b = histogram(v, normed=True, new=True)
        area = sum(a*diff(b))
        assert_almost_equal(area, 1)

        # Check with non constant bin width
        v = rand(n)*10
        bins = [0,1,5, 9, 10]
        a,b = histogram(v, bins, normed=True, new=True)
        area = sum(a*diff(b))
        assert_almost_equal(area, 1)


    def check_outliers_new(self):
        # Check that outliers are not tallied
        a = arange(10)+.5

        # Lower outliers
        h,b = histogram(a, range=[0,9], new=True)
        assert_equal(h.sum(),9)

        # Upper outliers
        h,b = histogram(a, range=[1,10], new=True)
        assert_equal(h.sum(),9)

        # Normalization
        h,b = histogram(a, range=[1,9], normed=True, new=True)
        assert_equal((h*diff(b)).sum(),1)

        # Weights
        w = arange(10)+.5
        h,b = histogram(a, range=[1,9], weights=w, normed=True, new=True)
        assert_equal((h*diff(b)).sum(),1)

        h,b = histogram(a, bins=8, range=[1,9], weights=w, new=True)
        assert_equal(h, w[1:-1])


    def check_type_new(self):
        # Check the type of the returned histogram
        a = arange(10)+.5
        h,b = histogram(a, new=True)
        assert(issubdtype(h.dtype, int))

        h,b = histogram(a, normed=True, new=True)
        assert(issubdtype(h.dtype, float))

        h,b = histogram(a, weights=ones(10, int), new=True)
        assert(issubdtype(h.dtype, int))

        h,b = histogram(a, weights=ones(10, float), new=True)
        assert(issubdtype(h.dtype, float))


    def check_weights_new(self):
        v = rand(100)
        w = ones(100)*5
        a,b = histogram(v,new=True)
        na,nb = histogram(v, normed=True, new=True)
        wa,wb = histogram(v, weights=w, new=True)
        nwa,nwb = histogram(v, weights=w, normed=True, new=True)
        assert_array_almost_equal(a*5, wa)
        assert_array_almost_equal(na, nwa)

        # Check weights are properly applied.
        v = linspace(0,10,10)
        w = concatenate((zeros(5), ones(5)))
        wa,wb = histogram(v, bins=arange(11),weights=w, new=True)
        assert_array_almost_equal(wa, w)

        # Check with integer weights
        wa, wb = histogram([1,2,2,4], bins=4, weights=[4,3,2,1], new=True)
        assert_array_equal(wa, [4,5,0,1])
        wa, wb = histogram([1,2,2,4], bins=4, weights=[4,3,2,1], normed=True, new=True)
        assert_array_equal(wa, array([4,5,0,1])/10./3.*4)

class TestHistogramdd(NumpyTestCase):
    def check_simple(self):
        x = array([[-.5, .5, 1.5], [-.5, 1.5, 2.5], [-.5, 2.5, .5], \
        [.5, .5, 1.5], [.5, 1.5, 2.5], [.5, 2.5, 2.5]])
        H, edges = histogramdd(x, (2,3,3), range = [[-1,1], [0,3], [0,3]])
        answer = asarray([[[0,1,0], [0,0,1], [1,0,0]], [[0,1,0], [0,0,1],
            [0,0,1]]])
        assert_array_equal(H,answer)
        # Check normalization
        ed = [[-2,0,2], [0,1,2,3], [0,1,2,3]]
        H, edges = histogramdd(x, bins = ed, normed = True)
        assert(all(H == answer/12.))
        # Check that H has the correct shape.
        H, edges = histogramdd(x, (2,3,4), range = [[-1,1], [0,3], [0,4]],
            normed=True)
        answer = asarray([[[0,1,0,0], [0,0,1,0], [1,0,0,0]], [[0,1,0,0],
            [0,0,1,0], [0,0,1,0]]])
        assert_array_almost_equal(H, answer/6., 4)
        # Check that a sequence of arrays is accepted and H has the correct
        # shape.
        z = [squeeze(y) for y in split(x,3,axis=1)]
        H, edges = histogramdd(z, bins=(4,3,2),range=[[-2,2], [0,3], [0,2]])
        answer = asarray([[[0,0],[0,0],[0,0]],
                          [[0,1], [0,0], [1,0]],
                          [[0,1], [0,0],[0,0]],
                          [[0,0],[0,0],[0,0]]])
        assert_array_equal(H, answer)

        Z = zeros((5,5,5))
        Z[range(5), range(5), range(5)] = 1.
        H,edges = histogramdd([arange(5), arange(5), arange(5)], 5)
        assert_array_equal(H, Z)

    def check_shape_3d(self):
        # All possible permutations for bins of different lengths in 3D.
        bins = ((5, 4, 6), (6, 4, 5), (5, 6, 4), (4, 6, 5), (6, 5, 4),
            (4, 5, 6))
        r = rand(10,3)
        for b in bins:
            H, edges = histogramdd(r, b)
            assert(H.shape == b)

    def check_shape_4d(self):
        # All possible permutations for bins of different lengths in 4D.
        bins = ((7, 4, 5, 6), (4, 5, 7, 6), (5, 6, 4, 7), (7, 6, 5, 4),
            (5, 7, 6, 4), (4, 6, 7, 5), (6, 5, 7, 4), (7, 5, 4, 6),
            (7, 4, 6, 5), (6, 4, 7, 5), (6, 7, 5, 4), (4, 6, 5, 7),
            (4, 7, 5, 6), (5, 4, 6, 7), (5, 7, 4, 6), (6, 7, 4, 5),
            (6, 5, 4, 7), (4, 7, 6, 5), (4, 5, 6, 7), (7, 6, 4, 5),
            (5, 4, 7, 6), (5, 6, 7, 4), (6, 4, 5, 7), (7, 5, 6, 4))

        r = rand(10,4)
        for b in bins:
            H, edges = histogramdd(r, b)
            assert(H.shape == b)

    def check_weights(self):
        v = rand(100,2)
        hist, edges = histogramdd(v)
        n_hist, edges = histogramdd(v, normed=True)
        w_hist, edges = histogramdd(v, weights=ones(100))
        assert_array_equal(w_hist, hist)
        w_hist, edges = histogramdd(v, weights=ones(100)*2, normed=True)
        assert_array_equal(w_hist, n_hist)
        w_hist, edges = histogramdd(v, weights=ones(100, int)*2)
        assert_array_equal(w_hist, 2*hist)

    def check_identical_samples(self):
        x = zeros((10,2),int)
        hist, edges = histogramdd(x, bins=2)
        assert_array_equal(edges[0],array([-0.5,  0. ,  0.5]))

class TestUnique(NumpyTestCase):
    def check_simple(self):
        x = array([4,3,2,1,1,2,3,4, 0])
        assert(all(unique(x) == [0,1,2,3,4]))
        assert(unique(array([1,1,1,1,1])) == array([1]))
        x = ['widget', 'ham', 'foo', 'bar', 'foo', 'ham']
        assert(all(unique(x) ==  ['bar', 'foo', 'ham', 'widget']))
        x = array([5+6j, 1+1j, 1+10j, 10, 5+6j])
        assert(all(unique(x) == [1+1j, 1+10j, 5+6j, 10]))

def compare_results(res,desired):
    for i in range(len(desired)):
        assert_array_equal(res[i],desired[i])

if __name__ == "__main__":
    NumpyTest().run()
