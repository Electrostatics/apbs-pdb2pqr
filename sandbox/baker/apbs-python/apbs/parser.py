""" APBS Python parser replacement """
from parameter import *
from input_file import *
import string
import re

# Define comment characters (can be multiple characters)
commentCharacters = "#"
# Define quote character (single character)
quoteCharacter = "\""

class Parser():
    """ Parse APBS input file """    
    def __init__(self):
        self.tokens = []
        self.inputFile = None

    def tokenize(self, data):
        """ Data must have readline() attribute """
        tokens = []
        for line in data:
            line = line.strip()
            # Get rid of comments
            regex = "[%s]+" % commentCharacters
            lineAndComments = re.split(regex, line)
            line = lineAndComments[0]
            # Split on quotes as if everything inside the quote was a token
            regex = "%s[^%s]+%s" % (quoteCharacter, quoteCharacter, quoteCharacter)
            while True:
                match = re.search(regex, line)
                if not match:
                    break
                start = match.start()
                end = match.end()
                tokens = tokens + line[:start].split()
                tokens = tokens + [line[(start+1):(end-1)]]
                line = line[end:]
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
    
