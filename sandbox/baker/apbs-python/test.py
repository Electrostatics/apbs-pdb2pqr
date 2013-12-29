""" APBS Python parser replacement """
from apbs import parser

if __name__ == '__main__':
    """ Test various aspects of the code """
    # Test mg-manual
    parser = parser.Parser()
    print "Testing mg-manual..."
    inpath = "examples/uber-input.in"
    infile = open(inpath, "rb")
    parser.feed(infile)
    inputFile = parser.parse()
    print inputFile
    
