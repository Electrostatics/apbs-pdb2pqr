""" APBS Python parser replacement """
from parameter import *
from input_file import *
import string

# Define comment character (single character)
commentCharacter = "#"

class Parser():
    """ Parse APBS input file """    
    def __init__(self):
        self.tokens = []
        self.inputFile = None

    def tokenize(self, data):
        """ Data must have readline() attribute """
        tokens = []
        while True:
            line = data.readline()
            if not line:
                break
            line = line.strip()
            lineAndComments = line.split(commentCharacter)
            line = lineAndComments[0]
            if line:
                tokens = tokens + line.split()
        return tokens
    
    def feed(self, data):
        """ Feed data into the parser. The data object must have a readline function. This function
        loads all of the data into the parser at once which is probably OK given how small APBS
        input files are. """
        self.tokens = self.tokens + self.tokenize(data)
 
    def parse(self):
        """ Parse the fed data, return InputFile object """
        self.inputFile = InputFile()
        self.inputFile.parse(self.tokens)
        return self.inputFile
    
