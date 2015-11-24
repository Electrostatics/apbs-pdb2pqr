""" Test APBS Python modules """
import unittest
import logging

logging.basicConfig(level=logging.DEBUG)

def print_test_results(result):
    """ Format and print results of tests """
    if result.wasSuccessful():
        print("Passed %s." % result)
    else:
        for error in result.errors:
            test_case = error[0]
            traceback = error[1]
            print(test_case, traceback)
        for failure in result.failures:
            print(failure)

def run_tests():
    """ Run tests """
    suite = unittest.defaultTestLoader.discover("apbs", "*.py")
    result = unittest.TestResult()
    suite.run(result)
    print_test_results(result)

if __name__ == "__main__":
    run_tests()
