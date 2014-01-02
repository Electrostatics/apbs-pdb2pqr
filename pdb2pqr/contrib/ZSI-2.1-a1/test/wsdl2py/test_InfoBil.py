#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import os, sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import FaultException
"""
Unittest 

WSDL:   http://javatest2.infodata.se/webservices/services/Infobil?wsdl
"""

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(TestCase, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(TestCase, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(TestCase, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(TestCase, 'test_'))
    return suite


class TestCase(ServiceTestCase):
    name = "test_InfoBil"
    client_file_name = "InfoBilServiceService_client.py"
    types_file_name  = "InfoBilServiceService_types.py"
    server_file_name = "InfoBilServiceService_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')

    def test_local_import(self):
        pass


if __name__ == "__main__" :
    main()

