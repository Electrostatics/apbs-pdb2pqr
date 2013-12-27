""" APBS input file class """
from read_section import *
from elec_section import *
from apolar_section import *
from print_section import *
from parameter import *


class InputFile(Parameter):
    """ Top-level APBS input file class """
    def __init__(self):
        self.tokens = None
        self.content_dict = { "read" : [], "elec" : [], "apolar" : [], "print" : []}
        
    @property
    def name(self):
        return "APBS INPUT FILE"
    
    def parse(self, tokens):
        """ This parses data read in with the feed() command """
        self.tokens = tokens
        token = self.tokens.pop(0)
        while True:
            sectionName = token.lower()
            if sectionName == "read":
                read = Read()
                read.parse(self.tokens)
                read.validate()
                self.content_dict["read"].append(read)
                print read
            elif sectionName == "elec":
                elec = Elec()
                elec.parse(self.tokens)
                elec.validate()
                self.content_dict["elec"].append(elec)
            else:
                errstr = "Unrecognized token %s" % token
                raise ValueError, errstr
            token = self.tokens.pop(0)
    
    def contents(self):
        return content_dict
    
    def validate(self):
        for values in self.content_dict.values():
            for param in values:
                param.validate()
    
    def __str__(self):
        outstr = ""
        for (key, value) in self.content_dict.items():
            outstr = outstr + "%s\n" % key
            outstr = outstr + "%s\n" % value
            outstr = outstr + "end\n"
        return outstr
        