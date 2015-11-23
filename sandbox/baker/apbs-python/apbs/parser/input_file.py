""" APBS input file class """
from .read_section import Read
from .print_section import Print
from .parameter import Parameter, Elec, Apolar

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
    def contents(self):
        return {"read" : self.reads, "calcs" : self.calcs, "print" : self.prints}
    def validate(self):
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
