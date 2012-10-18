#! /usr/bin/env python

from optparse import OptionParser
from datetime import datetime

apbs_binary_path = @APBS_BINARY_PATH@
inputgen_path = @INPUTGEN_PATH@

test_name = 'born'

input_dict = {
    'apbs-mol-auto'      : -2.297735411962E+02,
    'apbs-smol-auto'     : -2.290124171992E+02,
    'apbs-mol-parallel'  : -2.304918086635E+02, 
    'apbs-smol-parallel' : -2.293871354771E+02
     }
     
     
     
class Printer:
    
    def __init__( message_fd, logfile_fd ):
        self.message_fd = message_fd
        self.logfile_fd = logfile_fd
        
    def message( message ):
        print >> self.message_fd, message
        
    def log( message ):
        print >> self.logfile_fd, message
        


def process_serial( base_name ):
    
    input_name = '%s.in' % base_name
    output_name = '%s.out' % base_name
    output_file = open( output_name, 'w' )
    subprocess.call( "%s %s" % ( apbs_binary_path, input_name ), stdout = output_file )
    output_file = open( output_name, 'r' )
    output_text = output_file.read()
    output_pattern = r'Global net ELEC energy \= ([+-]?\d+\.\d+E[+-]\d+)'
    proc_result = float( re.search( output_pattern, output_text ).group( 1 ) )
    return proc_result
    


def process_parallel( base_name, procs, p ):

    p.message( "Splitting the input file into %d separate files using the inputgen utility" )
    p.message( "" )
    
    subprocess.call( "python %s --split %s" % ( inputgen_path, '%s.in' % base_name )
    result = 0.0
    
    for proc in range( procs ):
        proc_base_name = '%s-PE%d' % ( base_name, proc )
        proc_result = process_serial( proc_base_name, p )

        
        p.message( 'Processor %d result: %f' % ( proc, proc_result ) )
        p.message( '' )
        
        result += proc_result
        
    return result
        
        

def main():

    parser = OptionParser()
    parser.add_option(
        '-o', '--ocd', action='store_true', dest='ocd',
        help="Run APBS in OCD mode"
        )
    parser.add_option(
        '-l', '--log_file', dest='log_file', type='string', default='test.log',
        help="Save the test log to FILE.", metavar="FILE"
        )
    ( options, args ) = parser.parse_args()

    message_fd = sys.stdout
    logfile_fd = None
    try:
        logfile_fd = open( options.log_file, 'w' )
    except IOError as err:
        parser.error( "Could't open log_file %s: %s" % ( options.log_file, err.strerror ) )
    p = Printer( message_fd, logfile_fd )
    
    p.log( "Date:           %s" % str( datetime.now() ) )
    p.log( "Test Name:      %s" % test_name )
    p.log( "Target Results: %s" % ', '.join( [ str( r ) for r in input_dict.values() ] ) )
    
    
    
    for input_key in input_dict.keys():
        
        input_file = '%s.in' % input_key
        
        p.message( '----------------------------------------' )
        p.message( 'Testing input file %s' % input_file )
        p.message( '' )
        
        result = 0.0
        startime = datetime.now()
        match = re.search( r'\s*dime(\s+\d+)+', open( input_file, 'r' ).read() )
        if match != None:
            proc_list = [ int(p) for p in match.group( 1 ).split( ' ' ) ]
            procs = reduce( lambda x, y: x + y, proc_list )
            result = process_parallel( input_key, procs, p )
        else:
            result = process_serial( input_key )
            
        p.message( "Global net energy: %f" % result )
            
            
                
        
            
    
    
    
     
