""" APBS Python parser replacement. """

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

class CalcInput(parameter.Parameter):
    """ Top-level APBS input file class """
    def __init__(self):
        super(CalcInput, self).__init__()
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
    def parse(self, tokens):
        """ The input file has no parser """
        errstr = "This class should not be called for parsing."
        raise NotImplementedError(errstr)
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
    def encode(self, apbs_input):
        """ Encode an APBS CalcInput object into text """
        return apbs_input.__str__()

class TextDecoder:
    """ Parse APBS input file from plain text into Python object form """
    # Define comment characters (can be multiple characters)
    COMMENT_CHARACTERS = "#"
    # Define quote character (single character)
    QUOTE_CHARACTERS = "\""
    def __init__(self):
        self._tokens = []
        self.apbs_input = None
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
        """ Decode an APBS input file string and return an CalcInput object """
        with io.StringIO(input_string) as infile:
            self.feed(infile)
            return self.parse()
    def feed(self, data):
        """ Feed data into the parser. The data object must have a readline function. This function
        loads all of the data into the parser at once which is probably OK given how small APBS
        input files are. """
        self._tokens = self._tokens + self.tokenize(data)
    def parse(self):
        """ Parse the fed data (after self.feed() is called), return CalcInput object """
        apbs_input = CalcInput()
        token = self._tokens.pop(0)
        while True:
            section_name = token.lower()
            if section_name == "read":
                read = Read()
                read.parse(self._tokens)
                read.validate()
                apbs_input.reads.append(read)
            elif section_name == "elec":
                elec = Elec()
                elec.parse(self._tokens)
                elec.validate()
                apbs_input.calcs.append(elec)
            elif section_name == "apolar":
                apolar = Apolar()
                apolar.parse(self._tokens)
                apolar.validate()
                apbs_input.calcs.append(apolar)
            elif section_name == "print":
                print_object = Print()
                print_object.parse(self._tokens)
                print_object.validate()
                apbs_input.prints.append(print_object)
            elif section_name == "quit":
                break
            else:
                errstr = "Unrecognized token %s" % token
                raise ValueError(errstr)
            token = self._tokens.pop(0)
        return apbs_input

class JSONEncoder(json.JSONEncoder):
    """ Encode APBS CalcInput Python objects into JSON """
    def __init__(self, *args, **kwargs):
        super(JSONEncoder, self).__init__(*args, **kwargs)
    def default(self, obj):
        """ Encode object into JSON """
        if isinstance(obj, CalcInput):
            obj_list = []
            for read in obj.reads:
                obj_list.append({read.short_name() : read})
            for calc in obj.calcs:
                obj_list.append({calc.short_name() : calc})
            for print_ in obj.prints:
                obj_list.append({print_.short_name() : print_})
            return obj_list
        else:
            obj_dict = {}
            for key in obj.contents():
                value = getattr(obj, key)
                try:
                    obj_dict[value.short_name()] = value.get_value()
                except AttributeError:
                    obj_dict[key] = value
            return obj_dict

class JSONDecoder(json.JSONDecoder):
    """ Decode APBS JSON into CalcInput Python objects """
    def __init__(self, *args, **kwargs):
        super(JSONDecoder, self).__init__(*args, **kwargs)
        self.apbs_input = CalcInput()
    def decode(self, jstr, *args, **kwargs):
        jdict = super(JSONDecoder, self).decode(jstr)
        _LOGGER.debug(jdict)
        # TODO - I need to finish this JSON decoding class
        raise NotImplementedError

class _TestParser(unittest.TestCase):
    """ Test the parser """
    TEXT_INPUT_PATH = "./examples/uber-input.in"
    def setUp(self):
        self.text_input = open(self.TEXT_INPUT_PATH, "rt").read()
    def test_decode_text(self):
        """ Test text input """
        _LOGGER.info("Parsing TEXT input file...")
        TextDecoder().decode(self.text_input)
        _LOGGER.warn("This isn't a real unit test. Please review the output manually.")
    def test_encode_text(self):
        """ Test text output """
        apbs_from_text = TextDecoder().decode(self.text_input)
        _LOGGER.info("Writing TEXT output file...")
        text_output = TextEncoder().encode(apbs_from_text)
        _LOGGER.info("TEXT output:\n%s", text_output)
        _LOGGER.warn("This isn't a real unit test. Please review the output manually.")
    def test_encode_json(self):
        """ Test JSON output """
        _LOGGER.info("Writing JSON format output file...")
        apbs_from_text = TextDecoder().decode(self.text_input)
        json_input = json.dumps(apbs_from_text, indent=2, cls=JSONEncoder)
        _LOGGER.info("JSON output:\n%s", json_input)
        _LOGGER.warn("This isn't a real unit test. Please review the output manually.")
    @unittest.expectedFailure
    def test_decode_json(self):
        """ Test JSON input """
        # TODO - fix JSONDecoder so this test works
        apbs_from_text = TextDecoder().decode(self.text_input)
        json_input = json.dumps(apbs_from_text, indent=2, cls=JSONEncoder)
        _LOGGER.info("Parsing JSON input file...")
        _LOGGER.info("Result: %s", json.loads(json_input, cls=JSONDecoder))
        _LOGGER.warn("This isn't a real unit test. Please review the output manually.")
