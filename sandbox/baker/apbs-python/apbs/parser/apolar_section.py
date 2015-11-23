""" Parse the APOLAR input file section """
from .parameter import FLOAT_EPSILON, ParameterSection, OneStringParameter, OneFloatParameter
from . import elec_section

class Bconc(OneFloatParameter):
    """ This keyword specifies the bulk solvent density in &Aring;^-3. as described in the apolar
    calculation overview section. This coefficient multiplies the integral term of the apolar model
    discussed above and can be set to zero to eliminate integral contributions to the apolar
    solvation calculation. The syntax is

    bconc {density}

    where density is a floating point number giving the bulk solvent density in &Aring;^-3. """
    @property
    def name(self):
        return "bconc"
    def validate(self):
        if self.parm < 0:
            raise ValueError("bconc is negative")

class Name(OneStringParameter):
    """ The first (optional) argument is:

    name {id}

    where id is a unique string which can be assigned to the calculation to facilitate later
    operations (particularly in the PRINT statements). """
    @property
    def name(self):
        return "name"

class Calcenergy(elec_section.Calcenergy):
    """ This optional keyword controls energy output from an apolar solvation calculation. The
    syntax is

    calcenergy {flag}

    where flag is a string denoting what type of energy to calculate:

    * no - (Deprecated) don't calculate any energies.
    * total - Calculate and return total apolar energy for the entire molecule.
    * comps - Calculate and return total apolar energy for the entire molecule as well as the energy
    components for each atom.

    Note that this option must be used consistently for all calculations that will appear in
    subsequent PRINT statements. For example, if the statement

    print energy 1 - 2 end

    appears in the input file, then both calculations 1 and 2 must have calcenergy keywords present
    with the same values for flag. """
    pass

class Calcforce(elec_section.Calcforce):
    """ This optional keyword controls apolar force output. Its syntax is

    calcforce {flag}

    where flag is a string that specifies the types of force values to be returned:

    * no (Deprecated) - Don't calculate any forces
    * total - Calculate and return total apolar forces for the entire molecule
    * comps - Calculate and return total apolar forces for the entire molecule as well as force
    components for each atom.

    Note that this option must be used consistently for all calculations that will appear in
    subsequent PRINT statements. For example, if the statement

    print force 1 - 2 end

    appears in the input file, then both calculations 1 and 2 must have calcforce keywords present
    with the same values for flag. """
    pass

class Dpos(OneFloatParameter):
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
    @property
    def name(self):
        return "dpos"
    def validate(self):
        if self.parm < FLOAT_EPSILON:
            raise ValueError("dpos is zero or negative")

class Gamma(OneFloatParameter):
    """ This keyword specifies the surface tension coefficient for apolar solvation models.

    Its syntax is:

    gamma { value }

    where value is a floating point number designating the surface tension in units of kJ mol^-1
    &Aring;^-1. This term can be set to zero to eliminate SASA contributions to the apolar solvation
    calculations. """
    @property
    def name(self):
        return "gamma"
    def validate(self):
        if self.parm < 0:
            raise ValueError("negative gamma value")

class Grid(elec_section.Grid):
    """ This keyword specifies the quadrature grid spacing for volume integral calculations in
    apolar models.

    The syntax is:

    grid {hx hy hz}

    where hx hy hz are the quadrature spacings in the x-, y-, and z-directions in &Aring;. """
    pass

class Mol(elec_section.Mol):
    """ This term specifies the molecule for which the apolar calculation is to be performed.

    Its syntax is:

    mol {id}

    where id is the integer ID of the molecule for which the apolar calculation is to be performed.
    The molecule IDs are based on the order in which molecules are read by read mol statements,
    starting from 1. """
    pass

class Press(OneFloatParameter):
    """ This term specifies the solvent pressure p in kJ mol^-1 &Aring;^-3. This coefficient
    multiplies the volume term of the apolar model discussed here and can be set to zero to
    eliminate volume contributions to the apolar solvation calculation.

    The syntax is:

    press {value}

    where value is the floating point value of the pressure coefficient in kJ mol^-1 &Aring;-3. """
    @property
    def name(self):
        return "press"
    def validate(self):
        if self.parm < 0:
            raise ValueError("press is negative")

class Sdens(elec_section.Sdens):
    """ This keyword specifies the number of quadrature points per &Aring;^2 to use in surface terms
    (e.g., molecular surface, solvent accessible surface) for apolar calculations. The keyword is
    ignored when srad is 0.0 (e.g., for van der Waals surfaces) or when srfm is spl2 (e.g., for
    spline surfaces). The syntax is

    sdens {density}

    where density is a floating point number indicating the number of grid points per &Aring;^2.

    Users beware: there is a direct correlation between the value used for the sphere density, the
    accuracy of the results, and the APBS calculation time. """
    pass

class Srad(elec_section.Srad):
    """ This keyword specifies the radius of the solvent molecules; this parameter is used to define
    various solvent-related surfaces and volumes (see srfm). This value is usually set to 1.4
    &Aring; for a water-like molecular surface and set to 0 &Aring; for a van der Waals surface.

    Its syntax is

    srad {radius}

    where radius is the floating point value of the solvent radius (in &Aring;).

    This keyword is ignored for srfm spl2. """
    pass

class Srfm(OneStringParameter):
    """ This keyword specifies the model used to construct the solvent-related surface and volume.

    Its syntax is

    srfm {flag}

    where flag is a string that specifies the model used for surface and volume:

    * sacc - Solvent-accessible (also called "probe-inflated") surface and volume.

    Note that this keyword is under construction: we're in the process of adding additional surface
    definitions (e.g., spl2)."""
    def __init__(self):
        super(Srfm, self).__init__()
        self.allowed_values = ["sacc"]
    @property
    def name(self):
        return "srfm"

class Swin(elec_section.Swin):
    """ This keyword specifies the size of the support (i.e., the rate of change) for spline-based
    surface definitions (see srfm spl2). The value is usually set to 0.3 &Aring;.

    The syntax is

    swin {win}

    where win is the floating point value of the spline window (in &Aring;). """
    pass

class Temp(elec_section.Temp):
    """ This keyword specifies the temperature for the calculation. Its syntax is

    temp {T}

    where T is the floating point value of the temperature (in K) for the calculation. """
    pass

class Apolar(ParameterSection):
    """ APOLAR input file section

    This section is the main component for apolar solvation calculations in APBS runs. There may be
    several APOLAR sections, operating on different molecules or using different parameters for
    multiple runs on the same molecule. The syntax of this section is

    APOLAR [name id]
        {keywords...}
    END
    """
    def create_store_single_object(self, token_name, tokens):
        """ I'm not sure how safe this is but it sure saves a lot of code... This is designed to
        store parameters that only appear once per section block. """
        object_name = token_name.lower()
        class_name = object_name.capitalize()
        ctor = globals()[class_name]
        obj = ctor()
        obj.parse(tokens)
        print(class_name)
        obj.validate()
        self.content_dict[object_name] = obj
    def create_store_multiple_objects(self, token_name: str, tokens: list):
        """ I'm not sure how safe this is but it sure saves a lot of code... This is designed to
        store parameters that appear multiple times per section block. """
        object_name = token_name.lower()
        class_name = object_name.capitalize()
        ctor = globals()[class_name]
        obj = ctor()
        obj.parse(tokens)
        obj.validate()
        try:
            self.content_dict[object_name].append(obj)
        except KeyError:
            self.content_dict[object_name] = [obj]
    @property
    def name(self):
        return "apolar"
    def parse(self, tokens: list):
        token = tokens.pop(0)
        while True:
            token_name = token.lower()
            if token_name in ["name", "bconc", "calcenergy", "calcforce", "dpos", "gamma", "grid",
                              "mol", "press", "sdens", "srad", "srfm", "swin", "temp"]:
                # This section handles simple command statements
                self.create_store_single_object(token_name, tokens)
            elif token_name == "end":
                return
            else:
                errstr = "Unknown token (%s) for %s" % (token, self.name)
                raise ValueError(errstr)
            token = tokens.pop(0)
    def validate(self):
        required_values = ["name", "bconc", "dpos", "gamma", "grid",
                           "mol", "press", "sdens", "srad", "srfm", "swin", "temp"]
        for parm in required_values:
            if parm in self.content_dict:
                errstr = "APOLAR missing required parameter %s" % parm
                raise ValueError(errstr)
    def __str__(self) -> str:
        outstr = "apolar\n"
        keys = self.content_dict.keys()
        try:
            outstr = outstr + "\t%s\n" % self.content_dict["name"]
            keys.remove("name")
        except KeyError:
            pass
        for key in keys:
            values = self.content_dict[key]
            try:
                for value in values:
                    outstr = outstr + "\t%s\n" % value
            except TypeError:
                outstr = outstr + "\t%s\n" % values
        outstr = outstr + "end\n"
        return outstr
