""" Parse the PRINT input file section """
class Print:
    """ PRINT input file section """
    def __init__(self, tokens):
        self.tokens = tokens
        self.parse()
        
    def parse(self):
        """ Parse input file """
        for token in self.tokens:
            raise Exception, "Not implemented"
        