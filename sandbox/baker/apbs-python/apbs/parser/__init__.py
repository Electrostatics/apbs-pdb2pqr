""" APBS Python parser replacement. """
# TODO - eventually, this should be replaced with a JSON input format.

from . import parameter
from . import input_file
import string
import re
import unittest

import logging
_LOGGER = logging.getLogger("parser")

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
        self.input_file = input_file.InputFile()
        self.input_file.parse(self.tokens)
        return self.input_file

class _TestParser(unittest.TestCase):
    """ Test the parser """
    INPATH = "./examples/uber-input.in"
    def test_parser(self):
        """ Limited parser testing """
        self.parser = Parser()
        _LOGGER.info("Testing mg-manual...")
        infile = open(self.INPATH, "rt")
        self.parser.feed(infile)
        self.input_file = self.parser.parse()
        _LOGGER.info(self.input_file)
