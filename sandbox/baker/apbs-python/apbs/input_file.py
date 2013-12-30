""" APBS input file class """
from read_section import *
from elec_section import *
from apolar_section import *
from print_section import *
from parameter import *
import sys


class InputFile(Parameter):
    """ Top-level APBS input file class """
    def __init__(self):
        self.tokens = None
        self.reads = []
        self.calcs = []
        self.prints = []
        
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
                self.reads.append(read)
                print read
            elif sectionName == "elec":
                elec = Elec()
                elec.parse(self.tokens)
                elec.validate()
                self.calcs.append(elec)
                print elec
            elif sectionName == "apolar":
                apolar = Apolar()
                apolar.parse(self.tokens)
                apolar.validate()
                self.calcs.append(apolar)
            elif sectionName == "print":
                printObj = Print()
                printObj.parse(self.tokens)
                printObj.validate()
                self.prints.append(printObj)
            elif sectionName == "quit":
                return
            else:
                errstr = "Unrecognized token %s" % token
                raise ValueError, errstr
            token = self.tokens.pop(0)
    
    def contents(self):
        return { "read" : self.reads, "calcs" : self.calcs, "print" : self.prints }
    
    def validate(self):
        for values in self.content_dict.values():
            for param in values:
                param.validate()
    
    def __str__(self):
        outstr = ""
        for readSection in self.reads:
            outstr = outstr + "%s\n" % readSection
        for calcSection in self.calcs:
            outstr = outstr + "%s\n" % calcSection
        for printSection in self.prints:
            outstr = outstr + "%s\n" % printSection
        outstr = outstr + "quit\n"
        return outstr
        