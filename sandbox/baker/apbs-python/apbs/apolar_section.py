""" Parse the APOLAR input file section """
from parameter import *

class Name(OneStringParameter):
    """ The first (optional) argument is:

    name {id}
    
    where id is a unique string which can be assigned to the calculation to facilitate later
    operations (particularly in the PRINT statements). """
    @property
    def name(self):
        return "name"

class Apolar(ParameterSection):
    """ APOLAR input file section
    
    This section is the main component for apolar solvation calculations in APBS runs. There may be
    several APOLAR sections, operating on different molecules or using different parameters for
    multiple runs on the same molecule. The syntax of this section is
    
    APOLAR [name id]
        {keywords...}
    END
    """
    def createAndStoreSingleObject(self, tokenName, tokens):
        """ I'm not sure how safe this is but it sure saves a lot of code... This is designed to
        store parameters that only appear once per section block. """
        # TODO - The exec appears to execute in the current module namespace... can't figure out how
        # to change that to avoid duplicating code between classes.
        objectName = tokenName.lower()
        className = objectName.capitalize()
        execstr = "obj = %s()" % className
        exec execstr
        obj.parse(tokens)
        obj.validate()
        self.content_dict[objectName] = obj
    
    def createAndStoreMultipleObjects(self, tokenName, tokens):
        """ I'm not sure how safe this is but it sure saves a lot of code... This is designed to
        store parameters that appear multiple times per section block. """
        # TODO - The exec appears to execute in the current module namespace... can't figure out how
        # to change that to avoid duplicating code between classes.
        objectName = tokenName.lower()
        className = objectName.capitalize()
        execstr = "obj = %s()" % className
        exec execstr
        obj.parse(tokens)
        obj.validate()
        try:
            self.content_dict[objectName].append(obj)
        except KeyError:
            self.content_dict[objectName] = [obj]

    @property
    def name(self):
        return "apolar"
    
    def parse(self, tokens):
        token = tokens.pop(0)
        while True:
            tokenName = token.lower()
            if tokenName in ["name"]:
                # This section handles simple command statements
                self.createAndStoreSingleObject(tokenName, tokens)
            elif tokenName == "end":
                return
            else:
                errstr = "Unknown token (%s) for %s" % (token, self.name)
                raise ValueError, errstr
            token = tokens.pop(0)