""" APBS Python parser replacement. """
# TODO - eventually, this should be replaced with a JSON input format.

import re
import unittest
import logging
_LOGGER = logging.getLogger("parser")

from .elec_parser import Elec
from .apolar_parser import Apolar
from .print_parser import Print
from .read_parser import Read

class InputFile:
    """ Top-level APBS input file class """
    def __init__(self):
        self.tokens = None
        self.reads = []
        self.calcs = []
        self.prints = []
        self._short_name = "APBS INPUT FILE"
    def parse(self, tokens):
        """ This parses data read in with the feed() command """
        self.tokens = tokens
        token = self.tokens.pop(0)
        while True:
            section_name = token.lower()
            if section_name == "read":
                read = Read()
                read.parse(self.tokens)
                read.validate()
                self.reads.append(read)
            elif section_name == "elec":
                elec = Elec()
                elec.parse(self.tokens)
                elec.validate()
                self.calcs.append(elec)
            elif section_name == "apolar":
                apolar = Apolar()
                apolar.parse(self.tokens)
                apolar.validate()
                self.calcs.append(apolar)
            elif section_name == "print":
                print_object = Print()
                print_object.parse(self.tokens)
                print_object.validate()
                self.prints.append(print_object)
            elif section_name == "quit":
                return
            else:
                errstr = "Unrecognized token %s" % token
                raise ValueError(errstr)
            token = self.tokens.pop(0)
    def validate(self):
        """ Test the correctness of the object -- everything loaded correctly? """
        for values in self.reads + self.calcs + self.prints:
            for param in values:
                param.validate()
    def __str__(self):
        outstr = ""
        for read_section in self.reads:
            outstr = outstr + "%s\n" % read_section
        for calc_section in self.calcs:
            outstr = outstr + "%s\n" % calc_section
        for print_section in self.prints:
            outstr = outstr + "%s\n" % print_section
        outstr = outstr + "quit\n"
        return outstr

class Parser():
    """ Parse APBS input file """
    # Define comment characters (can be multiple characters)
    COMMENT_CHARACTERS = "#"
    # Define quote character (single character)
    QUOTE_CHARACTERS = "\""

    def __init__(self):
        self.tokens = []
        self.input_file = None

    def tokenize(self, data):
        """ Data must have readline() attribute """
        tokens = []
        for line in data:
            line = line.strip()
            # Get rid of comments
            regex = "[%s]+" % self.COMMENT_CHARACTERS
            line_comments = re.split(regex, line)
            line = line_comments[0]
            # Split on quotes as if everything inside the quote was a token
            regex = "%s[^%s]+%s" % (self.QUOTE_CHARACTERS, self.QUOTE_CHARACTERS,
                                    self.QUOTE_CHARACTERS)
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
        self.input_file = InputFile()
        self.input_file.parse(self.tokens)
        return self.input_file

class _TestParser(unittest.TestCase):
    """ Test the parser """
    INPATH = "./examples/uber-input.in"
    def __init__(self):
        super(_TestParser, self).__init__()
        self.parser = None
        self.input_file = None
    def test_parser(self):
        """ Limited parser testing """
        self.parser = Parser()
        _LOGGER.info("Testing mg-manual...")
        infile = open(self.INPATH, "rt")
        self.parser.feed(infile)
        self.input_file = self.parser.parse()
        _LOGGER.info(self.input_file)
