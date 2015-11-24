""" Handle the storage of APBS APOLAR block input file parameters """
import logging

_LOGGER = logging.getLogger("apolar-parser")

from . import parameter

class Name(parameter.Name):
    """ See base class documentation """
    pass
class Temp(parameter.Temp):
    """ See base class documentation """
    pass
class Calcenergy(parameter.Calcenergy):
    """ See base class documentation """
    pass
class Calcforce(parameter.Calcforce):
    """ See base class documentation """
    pass
class Mol(parameter.CalcMol):
    """ See base class documentation """
    pass
class Sdens(parameter.Sdens):
    """ See base class documentation """
    pass
class Srad(parameter.Srad):
    """ See base class documentation """
    pass
class Srfm(parameter.Srfm):
    """ See base class documentation """
    pass
class Swin(parameter.Swin):
    """ See base class documentation """
    pass
class Grid(parameter.Grid):
    """ See base class documentation """
    pass

class Bconc(parameter.OneFloatParameter):
    """ This keyword specifies the bulk solvent density in &Aring;^-3. as described in the apolar
    calculation overview section. This coefficient multiplies the integral term of the apolar model
    discussed above and can be set to zero to eliminate integral contributions to the apolar
    solvation calculation. The syntax is

    bconc {density}

    where density is a floating point number giving the bulk solvent density in &Aring;^-3. """
    def __init__(self):
        super(Bconc, self).__init__()
        self._short_name = "bconc"
    def validate(self):
        if self._parm < 0:
            raise ValueError("bconc is negative")

class Dpos(parameter.OneFloatParameter):
    """ This is the displacement used for finite-difference-based calculations of surface area
    derivatives.

    I know, this is a terrible way to calculate surface area derivatives -- we're working on
    replacing it with an analytic version. In the meantime, please use this parameter with caution.
    If anyone has code for a better method, please share!

    The syntax is

    dpos {displacement}

    where displacement is a floating point number indicating the finite difference displacement for
    force (surface area derivative) calculations in units of &Aring;.

    Important: this parameter is very dependent on sdens; e.g., smaller values of dpos require
    larger values of sdens. """
    def __init__(self):
        super(Dpos, self).__init__()
        self._short_name = "dpos"
    def validate(self):
        if self._parm < parameter.FLOAT_EPSILON:
            raise ValueError("dpos is zero or negative")

class Gamma(parameter.OneFloatParameter):
    """ This keyword specifies the surface tension coefficient for apolar solvation models.

    Its syntax is:

    gamma { value }

    where value is a floating point number designating the surface tension in units of kJ mol^-1
    &Aring;^-1. This term can be set to zero to eliminate SASA contributions to the apolar solvation
    calculations. """
    def __init__(self):
        super(Gamma, self).__init__()
        self._short_name = "gamma"
    def validate(self):
        if self._parm < 0:
            raise ValueError("negative gamma value")

class Press(parameter.OneFloatParameter):
    """ This term specifies the solvent pressure p in kJ mol^-1 &Aring;^-3. This coefficient
    multiplies the volume term of the apolar model discussed here and can be set to zero to
    eliminate volume contributions to the apolar solvation calculation.

    The syntax is:

    press {value}

    where value is the floating point value of the pressure coefficient in kJ mol^-1 &Aring;-3. """
    def __init__(self):
        super(Press, self).__init__()
        self._short_name = "press"
    def validate(self):
        if self._parm < 0:
            raise ValueError("press is negative")

class Apolar(parameter.ParameterSection):
    """ APOLAR input file section

    This section is the main component for apolar solvation calculations in APBS runs. There may be
    several APOLAR sections, operating on different molecules or using different parameters for
    multiple runs on the same molecule. The syntax of this section is

    APOLAR [name id]
        {keywords...}
    END
    """
    def __init__(self):
        super(Apolar, self).__init__()
        self._short_name = "apolar"
        self._allowed_keywords = {"name" : Name, "bconc" : Bconc, "calcenergy" : Calcenergy,
                                  "calcforce" : Calcforce, "dpos" : Dpos, "gamma" : Gamma,
                                  "grid" : Grid, "mol" : Mol, "press" : Press, "sdens" : Sdens,
                                  "srad" : Srad, "srfm" : Srfm, "swin" : Swin, "temp" : Temp}
    def parse(self, tokens: list):
        token = tokens.pop(0)
        while True:
            token_name = token.lower()
            if token_name in self._allowed_keywords:
                token_object = self._allowed_keywords[token_name]()
                token_object.parse(tokens)
                token_object.validate()
                setattr(self, token_name, token_object)
            elif token_name == "end":
                return
            else:
                errstr = "Unknown token (%s) for %s" % (token, self.short_name())
                raise ValueError(errstr)
            token = tokens.pop(0)
    def validate(self):
        required_values = {"name", "bconc", "dpos", "gamma", "grid",
                           "mol", "press", "sdens", "srad", "srfm", "swin", "temp"}
        for parm in required_values:
            if parm not in self.contents():
                errstr = "APOLAR missing required parameter %s" % parm
                _LOGGER.error("Incomplete parameters: %s", self.contents())
                raise ValueError(errstr)
    def __str__(self):
        outstr = self.format_block_start()
        keys = self.contents()
        keys.remove("name")
        for key in keys:
            values = getattr(self, key)
            try:
                for value in values:
                    outstr = outstr + "\t%s\n" % value
            except TypeError:
                outstr = outstr + "\t%s\n" % values
        outstr = outstr + "end\n"
        return outstr
