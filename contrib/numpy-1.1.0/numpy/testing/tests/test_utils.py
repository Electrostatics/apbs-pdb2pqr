import numpy as N
from numpy.testing.utils import *

import unittest


class _GenericTest(object):
    def _test_equal(self, a, b):
        self._assert_func(a, b)

    def _test_not_equal(self, a, b):
        passed = False
        try:
            self._assert_func(a, b)
            passed = True
        except AssertionError:
            pass

        if passed:
            raise AssertionError("a and b are found equal but are not")

    def test_array_rank1_eq(self):
        """Test two equal array of rank 1 are found equal."""
        a = N.array([1, 2])
        b = N.array([1, 2])

        self._test_equal(a, b)

    def test_array_rank1_noteq(self):
        """Test two different array of rank 1 are found not equal."""
        a = N.array([1, 2])
        b = N.array([2, 2])

        self._test_not_equal(a, b)

    def test_array_rank2_eq(self):
        """Test two equal array of rank 2 are found equal."""
        a = N.array([[1, 2], [3, 4]])
        b = N.array([[1, 2], [3, 4]])

        self._test_equal(a, b)

    def test_array_diffshape(self):
        """Test two arrays with different shapes are found not equal."""
        a = N.array([1, 2])
        b = N.array([[1, 2], [1, 2]])

        self._test_not_equal(a, b)

class TestEqual(_GenericTest, unittest.TestCase):
    def setUp(self):
        self._assert_func = assert_array_equal

    def test_generic_rank1(self):
        """Test rank 1 array for all dtypes."""
        def foo(t):
            a = N.empty(2, t)
            a.fill(1)
            b = a.copy()
            c = a.copy()
            c.fill(0)
            self._test_equal(a, b)
            self._test_not_equal(c, b)

        # Test numeric types and object
        for t in '?bhilqpBHILQPfdgFDG':
            foo(t)

        # Test strings
        for t in ['S1', 'U1']:
            foo(t)

    def test_generic_rank3(self):
        """Test rank 3 array for all dtypes."""
        def foo(t):
            a = N.empty((4, 2, 3), t)
            a.fill(1)
            b = a.copy()
            c = a.copy()
            c.fill(0)
            self._test_equal(a, b)
            self._test_not_equal(c, b)

        # Test numeric types and object
        for t in '?bhilqpBHILQPfdgFDG':
            foo(t)

        # Test strings
        for t in ['S1', 'U1']:
            foo(t)

    def test_nan_array(self):
        """Test arrays with nan values in them."""
        a = N.array([1, 2, N.nan])
        b = N.array([1, 2, N.nan])

        self._test_equal(a, b)

        c = N.array([1, 2, 3])
        self._test_not_equal(c, b)

    def test_string_arrays(self):
        """Test two arrays with different shapes are found not equal."""
        a = N.array(['floupi', 'floupa'])
        b = N.array(['floupi', 'floupa'])

        self._test_equal(a, b)

        c = N.array(['floupipi', 'floupa'])

        self._test_not_equal(c, b)

    def test_recarrays(self):
        """Test record arrays."""
        a = N.empty(2, [('floupi', N.float), ('floupa', N.float)])
        a['floupi'] = [1, 2]
        a['floupa'] = [1, 2]
        b = a.copy()

        self._test_equal(a, b)

        c = N.empty(2, [('floupipi', N.float), ('floupa', N.float)])
        c['floupipi'] = a['floupi'].copy()
        c['floupa'] = a['floupa'].copy()

        self._test_not_equal(c, b)


class TestAlmostEqual(_GenericTest, unittest.TestCase):
    def setUp(self):
        self._assert_func = assert_array_almost_equal


class TestRaises(unittest.TestCase):
    def setUp(self):
        class MyException(Exception):
            pass

        self.e = MyException

    def raises_exception(self, e):
        raise e

    def does_not_raise_exception(self):
        pass

    def test_correct_catch(self):
        f = raises(self.e)(self.raises_exception)(self.e)

    def test_wrong_exception(self):
        try:
            f = raises(self.e)(self.raises_exception)(RuntimeError)
        except RuntimeError:
            return
        else:
            raise AssertionError("should have caught RuntimeError")

    def test_catch_no_raise(self):
        try:
            f = raises(self.e)(self.does_not_raise_exception)()
        except AssertionError:
            return
        else:
            raise AssertionError("should have raised an AssertionError")

if __name__ == '__main__':
    unittest.main()
