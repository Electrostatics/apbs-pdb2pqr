import sys
from numpy.testing import *
set_package_path()
import numpy as np
restore_path()

def fft1(x):
    L = len(x)
    phase = -2j*np.pi*(np.arange(L)/float(L))
    phase = np.arange(L).reshape(-1,1) * phase
    return np.sum(x*np.exp(phase),axis=1)

class TestFFTShift(NumpyTestCase):
    def check_fft_n(self):
        self.failUnlessRaises(ValueError,np.fft.fft,[1,2,3],0)

class TestFFT1D(NumpyTestCase):
    def check_basic(self):
        rand = np.random.random
        x = rand(30) + 1j*rand(30)
        assert_array_almost_equal(fft1(x), np.fft.fft(x))

if __name__ == "__main__":
    NumpyTest().run()
