""" Handle the storage of APBS input file parameters """
from abc import ABCMeta, abstractmethod, abstractproperty

class Parameter:
    """ Parameter structure to simplify storing input file values """
    __metaclass__ = ABCMeta
    
    @abstractproperty
    def name(self):
        """ Each parameter class should define a name """
        pass
    
    @abstractmethod
    def parse(self, tokens):
        """ This function should pop tokens off the top of a stack and parse them.  The stack should
        be modified. """
        pass
    
    @abstractmethod
    def validate(self):
        """ Validate the contents of this parameter.  Raise Exception on error. """
        pass
    
    @abstractmethod
    def contents(self):
        """ Return the contents of this parameter as a dictionary """
        pass
    
    @abstractmethod
    def __str__(self):
        pass
    
class OneStringParameter(Parameter):
    """ Generic class for one-string parameter """
    def __init__(self):
        self.parm = None

    def parse(self, tokens):
        self.parm = tokens.pop(0)
    
    def validate(self):
        try:
            if not self.parm in self.allowed_values:
                errstr = "Unknown option %s for %s" % (self.parm, self.name)
                raise ValueError, errstr
        except AttributeError:
            pass
    
    def __str__(self):
        outstr = "%s %s\n" % (self.name, self.parm)
        return outstr
    
    def contents(self):
        return { self.name : self.parm }