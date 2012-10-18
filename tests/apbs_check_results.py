#! /usr/bin/env python

import sys, re
from apbs_logger import Logger

error_tolerance = 1e-4
from math import log10, floor

def round_sigfigs( x, sigfigs ):
    return round( x, sigfigs - int( floor( log10( abs( x ) ) ) ) -1 )

def check_results( computed_result, expected_result, input_file, logger, ocd ):
    if ocd:
        computed_result = round_sigfigs( computed_result, 12 )
        expected_result = round_sigfigs( expected_result, 12 )
    else:
        computed_result = round_sigfigs( computed_result, 7 )
        expected_result = round_sigfigs( expected_result, 7 )
    
    error = abs( ( computed_result - expected_result ) * 100.0 / 2.0 )
    
    if computed_result == expected_result:
        logger.message( "*** PASSED ***" )
        logger.log( "PASSED %.12e" % computed_result )
    elif error < error_tolerance:
        logger.message( "*** PASSED (with rounding error - see log) ***" )
        logger.log( "PASSED within error (%.12e; expected %.12e; %g%% error)" % ( computed_result, expected_result, error ) )
    else:
        logger.message( "*** FAILED ***" )
        logger.message( "   APBS returned      %.12e" % computed_result )
        logger.message( "   Expected result is %.12e (%g%% error)" % ( expected_result, error ) )
        logger.log( "FAILED (%.12e; expected %.12e; %g%% error)" % ( computed_result, expected_result, error ) )
    
    
    
if __name__ == '__main__':
    print >> sys.stderr, "The python source file %s is a module and not runnable" % sys.argv[ 0 ]
    sys.exit( 1 )
