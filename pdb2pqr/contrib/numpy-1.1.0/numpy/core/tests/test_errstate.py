# The following exec statement (or something like it) is needed to
# prevent SyntaxError on Python < 2.5. Even though this is a test,
# SyntaxErrors are not acceptable; on Debian systems, they block
# byte-compilation during install and thus cause the package to fail
# to install.

import sys
if sys.version_info[:2] >= (2, 5):
    exec """
from __future__ import with_statement
from numpy.core import *
from numpy.random import rand, randint
from numpy.testing import *



class TestErrstate(NumpyTestCase):


    def test_invalid(self):
        with errstate(all='raise', under='ignore'):
            a = -arange(3)
            # This should work
            with errstate(invalid='ignore'):
                sqrt(a)
            # While this should fail!
            try:
                sqrt(a)
            except FloatingPointError:
                pass
            else:
                self.fail()

    def test_divide(self):
        with errstate(all='raise', under='ignore'):
            a = -arange(3)
            # This should work
            with errstate(divide='ignore'):
                a / 0
            # While this should fail!
            try:
                a / 0
            except FloatingPointError:
                pass
            else:
                self.fail()

    def test_errcall(self):
        def foo(*args):
            print args
        olderrcall = geterrcall()
        with errstate(call=foo):
            assert(geterrcall() is foo), 'call is not foo'
            with errstate(call=None):
                assert(geterrcall() is None), 'call is not None'
        assert(geterrcall() is olderrcall), 'call is not olderrcall'

"""

if __name__ == '__main__':
    from numpy.testing import *
    NumpyTest().run()
