""" Parse the APOLAR input file section """
class Apolar:
    """ APOLAR input file section """
    def __init__(self, tokens):
        self.tokens = tokens
        self.parse()
        
    def parse(self):
        """ Parse input file """
        for token in self.tokens:
            raise Exception, "Not implemented"
        