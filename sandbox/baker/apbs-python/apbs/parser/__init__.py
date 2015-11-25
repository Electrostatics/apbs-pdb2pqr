""" APBS Python parser replacement. """
# TODO - eventually, this should be replaced with a JSON input format.

import json
import io
import re
import unittest
import logging
_LOGGER = logging.getLogger("parser")

from . import parameter
from .elec_parser import Elec
from .apolar_parser import Apolar
from .print_parser import Print
from .read_parser import Read

class APBSCalculation(parameter.Parameter):
    """ Top-level APBS input file class """
    def __init__(self):
        self._tokens = None
        self.reads = []
        self.calcs = []
        self.prints = []
        self._short_name = "APBS INPUT FILE"
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

class TextEncoder:
    """ Parse APBS input file objects into plain text """
    def encode(self, apbs_calculation):
        """ Encode an APBS APBSCalculation object into text """
        return apbs_calculation.__str__()

class TextDecoder:
    """ Parse APBS input file from plain text into Python object form """
    # Define comment characters (can be multiple characters)
    COMMENT_CHARACTERS = "#"
    # Define quote character (single character)
    QUOTE_CHARACTERS = "\""
    def __init__(self):
        self._tokens = []
        self.apbs_calculation = None
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
    def decode(self, input_string):
        """ Decode an APBS input file string and return an APBSCalculation object """
        with io.StringIO(input_string) as infile:
            self.feed(infile)
            return self.parse()
    def feed(self, data):
        """ Feed data into the parser. The data object must have a readline function. This function
        loads all of the data into the parser at once which is probably OK given how small APBS
        input files are. """
        self._tokens = self._tokens + self.tokenize(data)
    def parse(self):
        """ Parse the fed data (after self.feed() is called), return APBSCalculation object """
        apbs_calculation = APBSCalculation()
        token = self._tokens.pop(0)
        while True:
            section_name = token.lower()
            if section_name == "read":
                read = Read()
                read.parse(self._tokens)
                read.validate()
                apbs_calculation.reads.append(read)
            elif section_name == "elec":
                elec = Elec()
                elec.parse(self._tokens)
                elec.validate()
                apbs_calculation.calcs.append(elec)
            elif section_name == "apolar":
                apolar = Apolar()
                apolar.parse(self._tokens)
                apolar.validate()
                apbs_calculation.calcs.append(apolar)
            elif section_name == "print":
                print_object = Print()
                print_object.parse(self._tokens)
                print_object.validate()
                apbs_calculation.prints.append(print_object)
            elif section_name == "quit":
                break
            else:
                errstr = "Unrecognized token %s" % token
                raise ValueError(errstr)
            token = self._tokens.pop(0)
        return apbs_calculation

class JSONEncoder(json.JSONEncoder):
    """ Encode APBS APBSCalculation Python objects into JSON """
    def __init__(self, *args, **kwargs):
        super(JSONEncoder, self).__init__(*args, **kwargs)
    def default(self, obj):
        """ Encode object into JSON """
        obj_dict = {}
        for key in obj.contents():
            value = getattr(obj, key)
            try:
                obj_dict[value.short_name()] = value.get_value()
            except AttributeError:
                _LOGGER.warn("No specific handler for type %s", type(value))
                obj_dict[key] = value
        return obj_dict

class _TestParser(unittest.TestCase):
    """ Test the parser """
    INPATH = "./examples/uber-input.in"
    def setUp(self):
        self.input_string = open(self.INPATH, "rt").read()
        text_decoder = TextDecoder()
        _LOGGER.info("Parsing TEXT input file...")
        self.apbs_calculation = text_decoder.decode(self.input_string)
    def test_text(self):
        text_encoder = TextEncoder()
        _LOGGER.info("Input data in TEXT format:\n%s", text_encoder.encode(self.apbs_calculation))
        _LOGGER.warn("This isn't a real unit test. Please review the output manually.")
    def test_json(self):
        json_encoder = JSONEncoder()
        _LOGGER.info("Input data in JSON format:\n%s",
                     json.dumps(self.apbs_calculation, indent=2, cls=JSONEncoder))
