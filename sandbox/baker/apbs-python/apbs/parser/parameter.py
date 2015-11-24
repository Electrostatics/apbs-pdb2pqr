""" Handle the storage of APBS input file parameters """
# TODO - Replace globals()[class_name] with elegant decorator solution.
import sys
import logging

FLOAT_EPSILON = sys.float_info.epsilon
_LOGGER = logging.getLogger("parser")

# Set the following flag to True to use deprecated input format
USE_DEPRECATED = True

#### These are the generic classes for the parser

class Parameter(object):
    """ Parameter structure to simplify storing input file values """
    def __init__(self):
        self._short_name = None
    def short_name(self):
        """ Each parameter class should define a name """
        return self._short_name
    def parse(self, tokens):
        """ This function should pop tokens off the top of a stack and parse them.  The stack should
        be modified. """
        raise NotImplementedError
    def validate(self):
        """ Validate the contents of this parameter.  Raise Exception on error. """
        raise NotImplementedError
    def __str__(self):
        raise NotImplementedError
    def contents(self):
        """ Return a list of public variables for this class """
        var_list = []
        for key in vars(self).keys():
            if key[0] != "_":
                var_list.append(key)
        return var_list

class OneStringParameter(Parameter):
    """ Generic class for one-string parameter """
    def __init__(self):
        super(OneStringParameter, self).__init__()
        self._parm = None
        self._allowed_values = None
    def parse(self, tokens):
        self._parm = tokens.pop(0)
    def parm(self):
        return self._parm
    def validate(self):
        try:
            if not self._parm in self._allowed_values:
                errstr = "Unknown option %s for %s" % (self._parm, self.short_name())
                raise ValueError(errstr)
        except AttributeError:
            pass
    def __str__(self):
        outstr = "%s %s" % (self.short_name(), self._parm)
        return outstr

class OneIntegerParameter(Parameter):
    """ Generic class for one-integer parameter """
    def __init__(self):
        super(OneIntegerParameter, self).__init__()
        self._parm = None
    def parse(self, tokens):
        self._parm = int(tokens.pop(0))
    def parm(self):
        return self._parm
    def validate(self):
        if self._parm is None:
            errstr = "Missing value for parameter %s" % self.short_name()
            raise ValueError(errstr)
    def __str__(self):
        outstr = "%s %d" % (self.short_name(), self._parm)
        return outstr

class ParameterSection(Parameter):
    """ Complex parameters as found in an input file section """
    def __init__(self):
        super(ParameterSection, self).__init__()
    def format_block_start(self):
        """ Format the start of the block/section """
        outstr = self.short_name() + "\n"
        return outstr
    def create_store_multiple_objects(self, token_name, tokens):
        """ This is designed to store parameters that appear multiple times per section block. """
        object_name = token_name.lower()
        class_name = object_name.capitalize()
        ctor = globals()[class_name]
        obj = ctor()
        obj.parse(tokens)
        obj.validate()
        try:
            getattr(self, object_name).append(obj)
        except KeyError:
            setattr(self, object_name, [obj])
    def create_store_single_object(self, token_name, tokens):
        """ This is designed to store parameters that only appear once per section block. """
        object_name = token_name.lower()
        class_name = object_name.capitalize()
        ctor = globals()[class_name]
        obj = ctor()
        obj.parse(tokens)
        _LOGGER.debug("Creating object of class %s", class_name)
        obj.validate()
        setattr(self, object_name, obj)
    def parse(self, tokens):
        raise NotImplementedError
    def __str__(self):
        raise NotImplementedError

class OneFloatParameter(Parameter):
    """ Generic class for one-float parameter """
    def __init__(self):
        super(OneFloatParameter, self).__init__()
        self._parm = None
    def parse(self, tokens):
        self._parm = float(tokens.pop(0))
    def validate(self):
        if self._parm is None:
            errstr = "Missing value for parameter %s" % self.short_name
            raise ValueError(errstr)
    def __str__(self):
        outstr = "%s %g" % (self.short_name(), self._parm)
        return outstr

class ThreeFloatParameter(Parameter):
    """ Generic class for three-float parameter """
    def __init__(self):
        super(ThreeFloatParameter, self).__init__()
        self.xfloat = None
        self.yfloat = None
        self.zfloat = None
    def parse(self, tokens):
        self.xfloat = float(tokens.pop(0))
        self.yfloat = float(tokens.pop(0))
        self.zfloat = float(tokens.pop(0))
    def validate(self):
        # Validation happens through parsing
        pass
    def __str__(self):
        outstr = "%s %g %g %g" % (self.short_name(), self.xfloat, self.yfloat, self.zfloat)
        return outstr

class ThreeIntegerParameter(Parameter):
    """ Generic class for three-int parameter """
    def __init__(self):
        super(ThreeIntegerParameter, self).__init__()
        self.xint = None
        self.yint = None
        self.zint = None
    def parse(self, tokens):
        self.xint = int(tokens.pop(0))
        self.yint = int(tokens.pop(0))
        self.zint = int(tokens.pop(0))
    def validate(self):
        # Validation happens through parsing
        pass
    def __str__(self):
        outstr = "%s %d %d %d" % (self.short_name(), self.xint, self.yint, self.zint)
        return outstr

class FormatPathParameter(Parameter):
    """ Generic READ statement with format and path """
    def __init__(self):
        super(FormatPathParameter, self).__init__()
        self.format = None
        self.path = None
        self._allowed_values = None
    def parse(self, tokens):
        self.format = tokens.pop(0)
        path = "\"%s\"" % tokens.pop(0)
        self.path = path
    def validate(self):
        if not self.format in self._allowed_values:
            errstr = "Unknown token %s in %s" % (self.format, self.short_name)
            raise ValueError(errstr)
    def __str__(self):
        return " ".join([self.short_name(), self.format, self.path])

class Name(OneStringParameter):
    """ Usage: name {id}

    Since numerous ELEC/APOLAR blocks may appear in an APBS input file, it can be difficult to keep
    track of them all. It is possible to assign an optional name to each block to simplify the
    organizational process.

    {id} is an alphanumeric string denoting the "name" of the calculation block."""
    def __init__(self):
        super(Name, self).__init__()
        self._short_name = "name"
    def validate(self):
        if len(self._parm) == 0:
            raise ValueError("Can't have empty string parameter")

class Temp(OneFloatParameter):
    """ Specify the temperature for the Poisson-Boltzmann calculation. The syntax is:

    temp { T }

    where T is a floating point number indicating the temperature in K.

    Note that the temperature term is used for adjusting the ion distribution and scaling
    electrostatic potentials.  It is not used to model the temperature dependence of any
    dielectric terms. """
    def __init__(self):
        super(Temp, self).__init__()
        self._short_name = "temp"
    def validate(self):
        if self._parm < FLOAT_EPSILON:
            raise ValueError("temperature is zero or negative -- you violated the 3rd law!")

class Calcenergy(OneStringParameter):
    """ This optional keyword controls electrostatic energy output from a Poisson-Boltzmann
    calculation. Note that this option must be used consistently for all calculations that will
    appear in subsequent PRINT statements. For example, if the statement print energy 1 - 2 end
    appears in the input file, then both calculations 1 and 2 must have calcenergy keywords
    present with the same values for flag. The syntax for this keyword is:

    calcenergy { flag }

    where flag is a text string that specifies the types of energy values to be returned:

    * no - (Deprecated) don't calculate any energies.  This is the same as not including the
    calcenergy command in the input file.
    * total - Calculate and return total electrostatic energy for the entire molecule.
    * comps - Calculate and return total electrostatic energy for the entire molecule as well
    as electrostatic energy components for each atom. """
    def __init__(self):
        super(Calcenergy, self).__init__()
        self._allowed_values = ["no", "total", "comps"]
        self._short_name = "calcenergy"

class Calcforce(OneStringParameter):
    """ This optional keyword controls electrostatic force output from a multigrid (mg-manual,
    mg-auto, or mg-para) Poisson-Boltzmann calculation. Note that this option must be used
    consistently for all calculations that will appear in subsequent PRINT statements. For
    example, if the statement print force 1 - 2 end appears in the input file, then both
    calculations 1 and 2 must have calcforce keywords present with the same values for flag.
    The syntax for this keyword is:

    calcforce { flag }

    where flag is a text string that specifies the types of force values to be returned:

    * no - (Deprecated) don't calculate any forces.
    * total - Calculate and return total electrostatic and apolar forces for the entire
    molecule.
    * comps -     Calculate and return total electrostatic and apolar forces for the entire
    molecule as well as force components for each atom. """
    def __init__(self):
        super(Calcforce, self).__init__()
        self._allowed_values = ["no", "total", "comps"]
        self._short_name = "calcforce"

class Grid(ThreeFloatParameter):
    """ Specify the mesh grid spacings for multigrid mg-manual calculations.  This value may be
    different in each direction. The syntax is:

    grid {hx hy hz}

    where hx hy hz are the (floating point) grid spacings in the x-, y-, and z-directions
    (respectively) in &Aring;. """
    def __init__(self):
        super(Grid, self).__init__()
        self._short_name = "grid"
    def validate(self):
        if (self.xfloat < FLOAT_EPSILON) or (self.yfloat < FLOAT_EPSILON) \
        or (self.zfloat < FLOAT_EPSILON):
            errstr = "All grid spacings must be greater than zero!"
            raise ValueError(errstr)

class CalcMol(OneIntegerParameter):
    """ Specify the molecule for which the ELEC or APOLAR statement is to be evaluated. IDs are
    based on the order in which molecules are read by READ mol statements, starting from 1. The
    syntax is

    mol {id}

    where id is the integer ID of the molecule for which the Poisson-Boltzmann equation is to be
    solved. """
    def __init__(self):
        super(CalcMol, self).__init__()
        self._short_name = "mol"

class Sdens(OneFloatParameter):
    """ Specify the number of grid points per square-angstrom to use in discontinuous surface
    constructions (e.g., molecular surface and solvent-accessible surfaces). Ignored when srad is
    0.0 or srfm is spl2. There is a direct correlation between this value used for the surface
    sphere density, the accuracy of the surface calculations, and the APBS calculation time. The
    APBS "suggested" value is 10.0. The syntax of this command is

    sdens {density}

    where density is the floating point surface sphere density (in grid points/&Aring;^2).

    See also: srfm """
    def __init__(self):
        super(Sdens, self).__init__()
        self._short_name = "sdens"
    def validate(self):
        if self._parm < FLOAT_EPSILON:
            raise ValueError("sdens is zero or negative")

class Srad(OneFloatParameter):
    """ Specify the radius of the solvent molecules; this parameter is used to define the dielectric
    function for probe-based dielectric definitions (see srfm). This value is usually set to 1.4
    &Aring; for water. This keyword is ignored when any of the spline-based surfaces are used (e.g.,
    spl2, see srfm), since they are not probe-based. The syntax for this command is:

    srad {radius}

    where radius is the floating point solvent radius (in &Aring).

    See also: srfm """
    def __init__(self):
        super(Srad, self).__init__()
        self._short_name = "srad"
    def validate(self):
        if self._parm < FLOAT_EPSILON:
            return ValueError("srad is negative")

class Srfm(OneStringParameter):
    """ Specify the model used to construct the dielectric and ion-accessibility coefficients. The
    syntax is:

    srfm {flag}

    where flag is a string describing the coefficient model:

    * mol (ELEC only) - The dielectric coefficient is defined based on a molecular surface definition. The
    problem domain is divided into two spaces. The "free volume" space is defined by the union of
    solvent-sized spheres (see srad) which do not overlap with biomolecular atoms. This free volume
    is assigned bulk solvent dielectric values. The complement of this space is assigned
    biomolecular dielectric values. With a non-zero solvent radius (srad), this choice of
    coefficient corresponds to the traditional definition used for PB calculations. When the solvent
    radius is set to zero, this corresponds to a van der Waals surface definition. The ion-
    accessibility coefficient is defined by an "inflated" van der Waals model. Specifically, the
    radius of each biomolecular atom is increased by the radius of the ion species (as specified
    with the ion keyword). The problem domain is then divided into two spaces. The space inside the
    union of these inflated atomic spheres is assigned an ion-accessibility value of 0; the
    complement space is assigned bulk ion accessibility values.\
    * smol (ELEC only) - The dielectric and ion-accessibility coefficients are defined as for mol (see above).
    However, they are then "smoothed" by a 9-point harmonic averaging to somewhat reduce sensitivity
    to the grid setup as described by Bruccoleri et al. J Comput Chem 18 268-276, 1997 (journal web
    site).
    * spl2 (ELEC only) - The dielectric and ion-accessibility coefficients are defined by a cubic-spline surface
    as described by Im et al, Comp Phys Commun 111 (1-3) 59-75, 1998
    (doi:[10.1016/S0010-4655(98)00016-2). The width of the dielectric interface is controlled by the
    swin parameter.  These spline-based surface definitions are very stable with respect to grid
    parameters and therefore ideal for calculating forces. However, they require substantial
    reparameterization of the force field; interested users should consult Nina et al, Biophys Chem
    78 (1-2) 89-96, 1999 (doi:10.1016/S0301-4622(98)00236-1). Additionally, these surfaces can
    generate unphysical results with non-zero ionic strengths; this is an on-going area of
    development.
    * spl4 (ELEC only) - The dielectric and ion-accessibility coefficients are defined by a 7th order
    polynomial. This surface definition has characteristics similar to spl2, but provides higher
    order continuity necessary for stable force calculations with atomic multipole force fields (up
    to quadrupole).
    * sacc (APOLAR only) - Solvent-accessible (also called "probe-inflated") surface and volume. """
    # TODO - better to make this generic and then specialize derived classes for ELEC/APOLAR
    def __init__(self):
        super(Srfm, self).__init__()
        self._allowed_values = ["mol", "smol", "spl2", "spl4", "sacc"]
        self._short_name = "srfm"

class Swin(OneFloatParameter):
    """ Specify the size of the support (i.e., the rate of change) for spline-based surface
    definitions (see srfm). Usually 0.3 &Aring;. The syntax is:

    swin {win}

    where win is a floating point number for the spline window width (in &Aring;). Note that, per
    the analysis of Nina, Im, and Roux (doi:10.1016/S0301-4622(98)00236-1), the force field
    parameters (radii) generally need to be adjusted if the spline window is changed. """
    def __init__(self):
        super(Swin, self).__init__()
        self._short_name = "swin"
    def validate(self):
        if self._parm < FLOAT_EPSILON:
            raise ValueError("swin is zero or negative")
