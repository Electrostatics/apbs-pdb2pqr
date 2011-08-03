""" Test functions for limits module.
"""

from numpy.testing import *
import numpy.lib;reload(numpy.lib)
from numpy.lib.getlimits import finfo, iinfo
from numpy import single,double,longdouble
import numpy as np

##################################################

class TestPythonFloat(NumpyTestCase):
    def check_singleton(self):
        ftype = finfo(float)
        ftype2 = finfo(float)
        assert_equal(id(ftype),id(ftype2))

class TestSingle(NumpyTestCase):
    def check_singleton(self):
        ftype = finfo(single)
        ftype2 = finfo(single)
        assert_equal(id(ftype),id(ftype2))

class TestDouble(NumpyTestCase):
    def check_singleton(self):
        ftype = finfo(double)
        ftype2 = finfo(double)
        assert_equal(id(ftype),id(ftype2))

class TestLongdouble(NumpyTestCase):
    def check_singleton(self,level=2):
        ftype = finfo(longdouble)
        ftype2 = finfo(longdouble)
        assert_equal(id(ftype),id(ftype2))

class TestIinfo(NumpyTestCase):
    def check_basic(self):
        dts = zip(['i1', 'i2', 'i4', 'i8',
                   'u1', 'u2', 'u4', 'u8'],
                  [np.int8, np.int16, np.int32, np.int64,
                   np.uint8, np.uint16, np.uint32, np.uint64])
        for dt1, dt2 in dts:
            assert_equal(iinfo(dt1).min, iinfo(dt2).min)
            assert_equal(iinfo(dt1).max, iinfo(dt2).max)
        self.assertRaises(ValueError, iinfo, 'f4')

    def check_unsigned_max(self):
        types = np.sctypes['uint']
        for T in types:
            assert_equal(iinfo(T).max, T(-1))

if __name__ == "__main__":
    NumpyTest().run()
