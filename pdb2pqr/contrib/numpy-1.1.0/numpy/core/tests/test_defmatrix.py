from numpy.testing import *
set_package_path()
import numpy.core;reload(numpy.core)
from numpy.core import *
import numpy as np
restore_path()

class TestCtor(NumpyTestCase):
    def check_basic(self):
        A = array([[1,2],[3,4]])
        mA = matrix(A)
        assert all(mA.A == A)

        B = bmat("A,A;A,A")
        C = bmat([[A,A], [A,A]])
        D = array([[1,2,1,2],
                   [3,4,3,4],
                   [1,2,1,2],
                   [3,4,3,4]])
        assert all(B.A == D)
        assert all(C.A == D)

        vec = arange(5)
        mvec = matrix(vec)
        assert mvec.shape == (1,5)

class TestProperties(NumpyTestCase):
    def check_sum(self):
        """Test whether matrix.sum(axis=1) preserves orientation.
        Fails in NumPy <= 0.9.6.2127.
        """
        M = matrix([[1,2,0,0],
                   [3,4,0,0],
                   [1,2,1,2],
                   [3,4,3,4]])
        sum0 = matrix([8,12,4,6])
        sum1 = matrix([3,7,6,14]).T
        sumall = 30
        assert_array_equal(sum0, M.sum(axis=0))
        assert_array_equal(sum1, M.sum(axis=1))
        assert sumall == M.sum()

    def check_basic(self):
        import numpy.linalg as linalg

        A = array([[1., 2.],
                   [3., 4.]])
        mA = matrix(A)
        assert allclose(linalg.inv(A), mA.I)
        assert all(array(transpose(A) == mA.T))
        assert all(array(transpose(A) == mA.H))
        assert all(A == mA.A)

        B = A + 2j*A
        mB = matrix(B)
        assert allclose(linalg.inv(B), mB.I)
        assert all(array(transpose(B) == mB.T))
        assert all(array(conjugate(transpose(B)) == mB.H))

    def check_comparisons(self):
        A = arange(100).reshape(10,10)
        mA = matrix(A)
        mB = matrix(A) + 0.1
        assert all(mB == A+0.1)
        assert all(mB == matrix(A+0.1))
        assert not any(mB == matrix(A-0.1))
        assert all(mA < mB)
        assert all(mA <= mB)
        assert all(mA <= mA)
        assert not any(mA < mA)

        assert not any(mB < mA)
        assert all(mB >= mA)
        assert all(mB >= mB)
        assert not any(mB > mB)

        assert all(mA == mA)
        assert not any(mA == mB)
        assert all(mB != mA)

        assert not all(abs(mA) > 0)
        assert all(abs(mB > 0))

    def check_asmatrix(self):
        A = arange(100).reshape(10,10)
        mA = asmatrix(A)
        A[0,0] = -10
        assert A[0,0] == mA[0,0]

    def check_noaxis(self):
        A = matrix([[1,0],[0,1]])
        assert A.sum() == matrix(2)
        assert A.mean() == matrix(0.5)

class TestCasting(NumpyTestCase):
    def check_basic(self):
        A = arange(100).reshape(10,10)
        mA = matrix(A)

        mB = mA.copy()
        O = ones((10,10), float64) * 0.1
        mB = mB + O
        assert mB.dtype.type == float64
        assert all(mA != mB)
        assert all(mB == mA+0.1)

        mC = mA.copy()
        O = ones((10,10), complex128)
        mC = mC * O
        assert mC.dtype.type == complex128
        assert all(mA != mB)

class TestAlgebra(NumpyTestCase):
    def check_basic(self):
        import numpy.linalg as linalg

        A = array([[1., 2.],
                   [3., 4.]])
        mA = matrix(A)

        B = identity(2)
        for i in xrange(6):
            assert allclose((mA ** i).A, B)
            B = dot(B, A)

        Ainv = linalg.inv(A)
        B = identity(2)
        for i in xrange(6):
            assert allclose((mA ** -i).A, B)
            B = dot(B, Ainv)

        assert allclose((mA * mA).A, dot(A, A))
        assert allclose((mA + mA).A, (A + A))
        assert allclose((3*mA).A, (3*A))

class TestMatrixReturn(NumpyTestCase):
    def check_instance_methods(self):
        a = matrix([1.0], dtype='f8')
        methodargs = {
            'astype' : ('intc',),
            'clip' : (0.0, 1.0),
            'compress' : ([1],),
            'repeat' : (1,),
            'reshape' : (1,),
            'swapaxes' : (0,0)
            }
        excluded_methods = [
            'argmin', 'choose', 'dump', 'dumps', 'fill', 'getfield',
            'getA', 'getA1', 'item', 'nonzero', 'put', 'putmask', 'resize',
            'searchsorted', 'setflags', 'setfield', 'sort', 'take',
            'tofile', 'tolist', 'tostring', 'all', 'any', 'sum',
            'argmax', 'argmin', 'min', 'max', 'mean', 'var', 'ptp',
            'prod', 'std', 'ctypes', 'itemset'
            ]
        for attrib in dir(a):
            if attrib.startswith('_') or attrib in excluded_methods:
                continue
            f = eval('a.%s' % attrib)
            if callable(f):
                # reset contents of a
                a.astype('f8')
                a.fill(1.0)
                if attrib in methodargs:
                    args = methodargs[attrib]
                else:
                    args = ()
                b = f(*args)
                assert type(b) is matrix, "%s" % attrib
        assert type(a.real) is matrix
        assert type(a.imag) is matrix
        c,d = matrix([0.0]).nonzero()
        assert type(c) is matrix
        assert type(d) is matrix

class TestIndexing(NumpyTestCase):
    def check_basic(self):
        x = asmatrix(zeros((3,2),float))
        y = zeros((3,1),float)
        y[:,0] = [0.8,0.2,0.3]
        x[:,1] = y>0.5
        assert_equal(x, [[0,1],[0,0],[0,0]])

class TestNewScalarIndexing(NumpyTestCase):
    def setUp(self):
        self.a = matrix([[1, 2],[3,4]])

    def check_dimesions(self):
        a = self.a
        x = a[0]
        assert_equal(x.ndim, 2)

    def check_array_from_matrix_list(self):
        a = self.a
        x = array([a, a])
        assert_equal(x.shape, [2,2,2])

    def check_array_to_list(self):
        a = self.a
        assert_equal(a.tolist(),[[1, 2], [3, 4]])

    def check_fancy_indexing(self):
        a = self.a
        x = a[1, [0,1,0]]
        assert isinstance(x, matrix)
        assert_equal(x, matrix([[3,  4,  3]]))
        x = a[[1,0]]
        assert isinstance(x, matrix)
        assert_equal(x, matrix([[3,  4], [1, 2]]))
        x = a[[[1],[0]],[[1,0],[0,1]]]
        assert isinstance(x, matrix)
        assert_equal(x, matrix([[4,  3], [1,  2]]))

##     def check_vector_element(self):
##         x = matrix([[1,2,3],[4,5,6]])
##         assert_equal(x[0][0],1)
##         assert_equal(x[0].shape,(1,3))
##         assert_equal(x[:,0].shape,(2,1))

##         x = matrix(0)
##         assert_equal(x[0,0],0)
##         assert_equal(x[0],0)
##         assert_equal(x[:,0].shape,x.shape)

    def check_scalar_indexing(self):
        x = asmatrix(zeros((3,2),float))
        assert_equal(x[0,0],x[0][0])

    def check_row_column_indexing(self):
        x = asmatrix(np.eye(2))
        assert_array_equal(x[0,:],[[1,0]])
        assert_array_equal(x[1,:],[[0,1]])
        assert_array_equal(x[:,0],[[1],[0]])
        assert_array_equal(x[:,1],[[0],[1]])

    def check_boolean_indexing(self):
        A = arange(6)
        A.shape = (3,2)
        x = asmatrix(A)
        assert_array_equal(x[:,array([True,False])],x[:,0])
        assert_array_equal(x[array([True,False,False]),:],x[0,:])

    def check_list_indexing(self):
        A = arange(6)
        A.shape = (3,2)
        x = asmatrix(A)
        assert_array_equal(x[:,[1,0]],x[:,::-1])
        assert_array_equal(x[[2,1,0],:],x[::-1,:])



if __name__ == "__main__":
    NumpyTest().run()
