""" APBS Python parser replacement. """

import json
import io
import re
import unittest
import logging
_LOGGER = logging.getLogger(__name__)

from . import parameter
from .elec_parser import Elec
from .apolar_parser import Apolar
from .print_parser import Print
from .read_parser import Read

class CalcInput(parameter.Parameter):
    """ Top-level APBS input file class """
    def __init__(self, *args, **kwargs):
        super(CalcInput, self).__init__(*args, **kwargs)
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
            self.apbs_input = apbs_input
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
    read_texts = [["read", "charge gz charge.dx.gz", "diel gz dielx.dx.gz diely.dx.gz dielz.dx.gz",
                  "kappa gz \"./This path has spaces/kappa.dx.gz\"", "mesh mcsf mesh.m",
                  "mol pqr mol.pqr", "pot gz pot.dx.gz", "end"]]
    calc_texts = [["elec", "name mg-auto1", "mg-auto", "bcfl mdh # Make sure inline comments work...",
                   "calcforce comps", "chgm spl0", "cgcent mol 1", "cglen 100.1 91.6 55.4",
                   "dime 65 65 65", "fgcent 0.1 -0.009 10.0", "fglen 10.0 9.1 5.5",
                   "ion charge +1 conc 0.1 radius 2.0", "ion charge -1 conc 0.1 radius 2.0", "lpbe",
                   "mol 1", "pdie 2.0", "sdens 10.0", "sdie 78.54", "srad 1.4", "srfm smol",
                   "swin 0.2", "temp 298.15", "write pot gz pot", "writemat poisson mat", "end"],
                  ["elec", "name mg-para2", "mg-para", "async 4", "bcfl mdh", "calcenergy total",
                   "calcforce comps", "chgm spl0", "cgcent +1.4 -0.133 0.1",
                   "cglen 110.0 95.0 62.3", "dime 97 129 130", "etol 0.1", "fgcent mol 1",
                   "fglen 11.0 9.5 6.2", "ion charge +3 conc 0.2 radius 3.0",
                   "ion charge -3 conc 0.2 radius 3.0", "npbe", "mol 1", "ofrac 0.2", "pdie 2.0",
                   "pdime 5 5 5", "sdens 20.0", "sdie 80.0", "srad 1.4", "srfm spl2", "swin 0.3",
                   "temp 298.15", "usemap kappa 1", "write kappa dx \"kappa.dx.gz\"", "end"],
                  ["elec", "name mg-manual3", "mg-manual", "bcfl sdh", "calcenergy no",
                   "calcforce no", "chgm spl4", "dime 161 97 129", "etol 0.1", "gcent mol 1",
                   "grid 0.51 0.52 0.53", "lpbe", "mol 1", "nlev 3", "pdie 2.0", "sdens 30.0",
                   "sdie 80.0", "srad 1.4", "srfm mol", "swin 0.2", "temp 50.0",
                   "write pot gz \"./directory with spaces/pot\"", "end"],
                  ["elec", "name fe-manual", "fe-manual", "akeyPRE unif", "akeySOLVE resi",
                   "bcfl zero", "calcenergy comps", "calcforce total", "chgm spl2",
                   "domainLength 50.0 39.2 55.2", "ekey simp", "etol 0.4",
                   "ion charge +2 conc 0.25 radius 1.0", "ion charge -1 conc 0.50 radius 1.0",
                   "nrpbe", "maxsolve 100", "maxvert 10000", "mol 1", "pdie 2.0", "sdens 50.0",
                   "sdie 80.0", "srad 1.4", "srfm smol", "swin 0.3", "temp 300.0", "end"],
                  ["elec", "name mg-dummy", "mg-manual", "bcfl sdh", "calcenergy no",
                   "calcforce no", "chgm spl4", "dime 161 97 129", "etol 0.1", "gcent mol 1",
                   "grid 0.511 0.522 0.533", "lpbe", "mol 1", "nlev 3", "pdie 2.0", "sdens 30.0",
                   "sdie 80.0", "srad 1.4", "srfm mol", "swin 0.2", "temp 50.0",
                   "write pot gz \"./directory with spaces/pot\"", "end"],
                  ["apolar", "name apolar", "bconc 0.1", "calcenergy total", "calcforce total",
                   "dpos 0.1", "gamma 25.0 # Test inline comment", "grid 0.11 0.12 0.13", "mol 1",
                   "press 0.2", "sdens 20.1", "srad 1.2", "srfm sacc", "swin 0.3", "temp 298",
                   "end"]]
    print_texts = [["print", "elecEnergy 1 - 2", "end"]]
    def setUp(self):
        self.text_input = ""
        for read in self.read_texts:
            self.text_input = self.text_input + "\n" + "\n".join(read)
        for calc in self.calc_texts:
            self.text_input = self.text_input + "\n" + "\n".join(calc)
        for print_ in self.print_texts:
            self.text_input = self.text_input + "\n" + "\n".join(print_)
        self.text_input = self.text_input + "\nquit\n"
        _LOGGER.debug("*** START TEST INPUT ***\n%s\n *** END TEST INPUT ***", self.text_input)
    def compare_results(self, start_indices, end_indices, text_output_lines, reference_lists):
        import difflib
        junk_chars = lambda x : x in " \""
        self.assertTrue(junk_chars(" "))
        self.assertTrue(junk_chars("\""))
        self.assertFalse(junk_chars("a"))
        matcher = difflib.SequenceMatcher(isjunk=junk_chars)
        for (isection, start_index) in enumerate(start_indices):
            end_index = end_indices.pop(0)
            test_section = sorted(text_output_lines[start_index:end_index+1])
            ref_section = sorted(reference_lists[isection])
            test_section_lines = {}
            for line in test_section:
                keyword = line.split()[0].lower()
                # Not testing for proper quotation behavior because of variation in reference input
                line = line.replace("\"", "")
                try:
                    test_section_lines[keyword].append(line)
                except KeyError:
                    test_section_lines[keyword] = [line]
            ref_section_lines = {}
            for line in ref_section:
                keyword = line.split()[0].lower()
                line = line.replace("\"", "")
                # Remove comments
                try:
                    line = line[:line.index("#")].strip()
                except ValueError:
                    pass
                try:
                    ref_section_lines[keyword].append(line)
                except KeyError:
                    ref_section_lines[keyword] = [line]
            self.assertTrue(len(test_section) == len(ref_section))
            for keyword in ref_section_lines.keys():
                if keyword in ["lpbe", "npbe", "nrpbe", "lrpbe"]:
                    test_lines = sorted(test_section_lines["eqntype"])
                    ref_lines = sorted(["eqntype %s" % line for line in ref_section_lines[keyword]])
                elif keyword in ["mg-auto", "mg-para", "mg-manual", "fe-manual"]:
                    test_lines = sorted(test_section_lines["solvtype"])
                    ref_lines = sorted(["solvtype %s" % line for line in ref_section_lines[keyword]])
                elif keyword == "ion":
                    normalize_charge = lambda line : float(line.split()[2])
                    test_lines = sorted(test_section_lines[keyword], key=normalize_charge)
                    ref_lines = sorted(ref_section_lines[keyword], key=normalize_charge)
                else:
                    test_lines = sorted(test_section_lines[keyword])
                    ref_lines = sorted(ref_section_lines[keyword])
                matcher.set_seqs(test_lines, ref_lines)
                for opcode in matcher.get_opcodes():
                    if opcode[0] != "equal":
                        # Most errors are due to rounding/printf ambiguity
                        if keyword in ["pdie", "sdie", "temp", "sdens", "gamma"]:
                            test_val = float(test_lines[0].split()[1])
                            ref_val = float(ref_lines[0].split()[1])
                            self.assertTrue(test_val == ref_val)
                        elif keyword in ["glen", "fglen", "cglen", "cgcent", "fgcent", "gcent",
                                         "domainlength"]:
                            for idim in range(3):
                                test_val = float(test_lines[0].split()[idim+1])
                                ref_val = float(ref_lines[0].split()[idim+1])
                                self.assertTrue(test_val == ref_val)
                        elif keyword == "ion":
                            for iline, ref_line in enumerate(ref_lines):
                                test_line = test_lines[iline]
                                for idim in [2, 4, 6]:
                                    test_val = float(test_line.split()[idim])
                                    ref_val = float(ref_line.split()[idim])
                                    self.assertTrue(test_val == ref_val)
                        elif keyword == "dime":
                            # This tests that the input dime with strange values were correctly
                            # adjusted by the input parsing code
                            test_line = test_lines[0]
                            self.assertTrue(test_line == "dime 97 129 87")
                            ref_line = ref_lines[0]
                            self.assertTrue(ref_line == "dime 97 129 130")
                        elif keyword == "domainLength":
                            test_line = test_lines[0].lower()
                            ref_line = ref_lines[0].lower()
                            self.assertTrue(test_line == ref_line)
                        else:
                            _LOGGER.debug(opcode)
                            _LOGGER.debug(keyword)
                            _LOGGER.debug(test_lines, ref_lines)
                            _LOGGER.debug(matcher.find_longest_match(0, len(test_lines), 0, len(ref_lines)))
                            self.assertTrue(test_lines[0].lower() == ref_lines[0].lower())

    def test_valid_text(self):
        """ Test valid text input """
        # This part loads the text
        _LOGGER.info("Testing APBS text input file parsing...")
        calc_input = TextDecoder().decode(self.text_input)
        # This creates the output for comparison with input
        _LOGGER.info("Creating TEXT output file...")
        text_output = TextEncoder().encode(calc_input)
        # This compares the input and output
        text_output_lines = []
        import string
        whitespace = string.whitespace + "\n\r"
        for line in io.StringIO(text_output).readlines():
            line = line.strip(whitespace)
            if len(line) > 0:
                text_output_lines.append(line)
        read_indices = [n for (n, e) in enumerate(text_output_lines) if e == "read"]
        calc_indices = [n for (n, e) in enumerate(text_output_lines) if e in ["elec", "apolar"]]
        print_indices = [n for (n, e) in enumerate(text_output_lines) if e == "print"]
        end_indices = [n for (n, e) in enumerate(text_output_lines) if e == "end"]
        # Compare the read sections
        _LOGGER.info("Testing READ parsing...")
        self.compare_results(read_indices, end_indices, text_output_lines, self.read_texts)
        _LOGGER.info("Testing ELEC/APOLAR parsing...")
        self.compare_results(calc_indices, end_indices, text_output_lines, self.calc_texts)
    def test_encode_json(self):
        """ Test JSON output """
        _LOGGER.info("Testing APBS JSON output file generation...")
        _LOGGER.warn("This isn't a real unit test -- sorry.")
        apbs_from_text = TextDecoder().decode(self.text_input)
        json_input = json.dumps(apbs_from_text, indent=2, cls=JSONEncoder)
        _LOGGER.debug("JSON output:\n%s", json_input)
