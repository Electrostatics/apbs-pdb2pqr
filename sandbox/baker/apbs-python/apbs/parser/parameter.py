""" Handle the storage of APBS input file parameters """
# TODO - Replace globals()[class_name] with elegant decorator solution.
# TODO - Break this back into multiple source files -- it's too long.
import sys
import logging

from .utility import factors, product

_LOGGER = logging.getLogger("parser")
FLOAT_EPSILON = sys.float_info.epsilon

#### These are the generic classes

class Parameter(object):
    """ Parameter structure to simplify storing input file values """
    def name(self):
        """ Each parameter class should define a name """
        raise NotImplementedError
    def parse(self, tokens):
        """ This function should pop tokens off the top of a stack and parse them.  The stack should
        be modified. """
        raise NotImplementedError
    def validate(self):
        """ Validate the contents of this parameter.  Raise Exception on error. """
        raise NotImplementedError
    def contents(self):
        """ Return the contents of this parameter as a dictionary """
        raise NotImplementedError
    def __str__(self):
        raise NotImplementedError

class OneStringParameter(Parameter):
    """ Generic class for one-string parameter """
    def __init__(self):
        self.parm = None
        self.allowed_values = None
    def name(self):
        raise NotImplementedError
    def parse(self, tokens):
        self.parm = tokens.pop(0)
    def validate(self):
        try:
            if not self.parm in self.allowed_values:
                errstr = "Unknown option %s for %s" % (self.parm, self.name())
                raise ValueError(errstr)
        except AttributeError:
            pass
    def __str__(self):
        outstr = "%s %s" % (self.name, self.parm)
        return outstr
    def contents(self):
        return {self.name : self.parm}

class OneIntegerParameter(Parameter):
    """ Generic class for one-integer parameter """
    def __init__(self):
        self.parm = None
    def name(self):
        raise NotImplementedError
    def parse(self, tokens):
        self.parm = int(tokens.pop(0))
    def validate(self):
        if self.parm is None:
            errstr = "Missing value for parameter %s" % self.name
            raise ValueError(errstr)
    def __str__(self):
        outstr = "%s %d" % (self.name, self.parm)
        return outstr
    def contents(self):
        return {self.name : self.parm}

class ParameterSection(Parameter):
    """ Complex parameters as found in an input file section """
    def __init__(self):
        self.content_dict = {}
        self.my_name = None
    def name(self):
        if not self.my_name:
            raise NotImplementedError
        else:
            return self.my_name
    def format_block_start(self):
        """ Format the start of the block/section """
        outstr = self.name() + "\n"
        keys = list(self.content_dict.keys())
        try:
            outstr = outstr + "\t%s\n" % self.content_dict["name"]
            keys.remove("name")
        except KeyError:
            pass
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
            self.content_dict[object_name].append(obj)
        except KeyError:
            self.content_dict[object_name] = [obj]
    def create_store_single_object(self, token_name, tokens):
        """ This is designed to store parameters that only appear once per section block. """
        object_name = token_name.lower()
        class_name = object_name.capitalize()
        ctor = globals()[class_name]
        obj = ctor()
        obj.parse(tokens)
        _LOGGER.debug("Creating object of class %s", class_name)
        obj.validate()
        self.content_dict[object_name] = obj
    def parse(self, tokens):
        raise NotImplementedError
    def __str__(self):
        raise NotImplementedError
    def contents(self):
        return self.content_dict

class OneFloatParameter(Parameter):
    """ Generic class for one-float parameter """
    def __init__(self):
        self.parm = None
    def name(self):
        raise NotImplementedError
    def parse(self, tokens):
        self.parm = float(tokens.pop(0))
    def validate(self):
        if self.parm is None:
            errstr = "Missing value for parameter %s" % self.name
            raise ValueError(errstr)
    def __str__(self):
        outstr = "%s %g" % (self.name, self.parm)
        return outstr
    def contents(self):
        return {self.name : self.parm}

class ThreeFloatParameter(Parameter):
    """ Generic class for three-float parameter """
    def __init__(self):
        self.xfloat = None
        self.yfloat = None
        self.zfloat = None
    def name(self):
        raise NotImplementedError
    def parse(self, tokens):
        self.xfloat = float(tokens.pop(0))
        self.yfloat = float(tokens.pop(0))
        self.zfloat = float(tokens.pop(0))
    def validate(self):
        # Validation happens through parsing
        pass
    def __str__(self):
        outstr = "%s %g %g %g" % (self.name, self.xfloat, self.yfloat, self.zfloat)
        return outstr
    def contents(self):
        return {"xfloat" : self.xfloat, "yfloat" : self.yfloat, "zfloat" : self.zfloat}

class ThreeIntegerParameter(Parameter):
    """ Generic class for three-int parameter """
    def __init__(self):
        self.xint = None
        self.yint = None
        self.zint = None
    def name(self):
        raise NotImplementedError
    def parse(self, tokens):
        self.xint = int(tokens.pop(0))
        self.yint = int(tokens.pop(0))
        self.zint = int(tokens.pop(0))
    def validate(self):
        # Validation happens through parsing
        pass
    def __str__(self):
        outstr = "%s %d %d %d" % (self.name, self.xint, self.yint, self.zint)
        return outstr
    def contents(self):
        return {"xint" : self.xint, "yint" : self.yint, "zint" : self.zint}

class FormatPathParameter(Parameter):
    """ Generic READ statement with format and path """
    def __init__(self):
        super(FormatPathParameter, self).__init__()
        self.format = None
        self.path = None
        self.allowed_values = None
    def name(self):
        raise NotImplementedError
    def contents(self):
        return {"format" : self.format, "path" : self.path}
    def parse(self, tokens):
        self.format = tokens.pop(0)
        path = "\"%s\"" % tokens.pop(0)
        self.path = path
    def validate(self):
        if not self.format in self.allowed_values:
            errstr = "Unknown token %s in %s" % (self.format, self.name)
            raise ValueError(errstr)
    def __str__(self):
        return " ".join([self.name, self.format, self.path])

class Name(OneStringParameter):
    """ Usage: name {id}

    Since numerous ELEC/APOLAR blocks may appear in an APBS input file, it can be difficult to keep
    track of them all. It is possible to assign an optional name to each block to simplify the
    organizational process.

    {id} is an alphanumeric string denoting the "name" of the calculation block."""
    def __init__(self):
        super(Name, self).__init__()
        self.my_name = "name"
    def validate(self):
        if len(self.parm) == 0:
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
        self.my_name = "temp"
    def validate(self):
        if self.parm < FLOAT_EPSILON:
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
        self.allowed_values = ["no", "total", "comps"]
    @property
    def name(self):
        return "calcenergy"

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
        self.allowed_values = ["no", "total", "comps"]
    @property
    def name(self):
        return "calcforce"

### These are the ELEC-specific classes

class Glen(ThreeFloatParameter):
    """ Specify the mesh domain lengths for multigrid mg-manual calculations.  These lengths may be
    different in each direction. The syntax is:

    glen {xlen ylen zlen}

    where xlen ylen zlen are the (floating point) grid lengths in the x-, y-, and z-directions
    (respectively) in &Aring;.

    See also: grid  """
    def __init__(self):
        super(Glen, self).__init__()
    @property
    def name(self):
        return "glen"
    def validate(self):
        if (self.xfloat < FLOAT_EPSILON) or (self.yfloat < FLOAT_EPSILON) or (self.zfloat <
                                                                              FLOAT_EPSILON):
            errstr = "One of the grid lengths is zero or negative (%g, %g, %g)" % (self.xfloat,
                                                                                   self.yfloat,
                                                                                   self.zfloat)
            raise ValueError(errstr)

class Solvtype(OneStringParameter):
    """ Usage:  solvtype {type}

    Specify the type of electrostatics calculation to be performed where {type} is one of:

    * fe-manual - manual finite element method
    * mg-auto - automatic multigrid setup
    * mg-dummy - dummy multigrid setup for writing out maps
    * mg-manual - manual multigrid setup
    * mg-para - parallel multigrid setup

    This isn't actually official yet.  I would like to change the default behavior of APBS
    from stating the ELEC calculation type without any parameters to adding the "type" keyword."""
    def __init__(self):
        super(Solvtype, self).__init__()
        self.allowed_values = ["mg-auto", "mg-para", "mg-manual", "mg-dummy", "fe-manual"]
    @property
    def name(self):
        return "solvtype"

class Akeypre(OneStringParameter):
    """ Usage: akeyPRE {key}

    Specifies how the initial finite element mesh should be constructed (from refinement of a
    very coarse 8-tetrahedron mesh prior to the solve-estimate-refine iteration in fe-manual
    finite element calculations.

    key is a text string that specifies the method used to guide initial refinement and takes
    one of the values:

    * unif - Uniform refinement
    * geom - Geometry-based refinement at molecular surfaces and charges """
    def __init__(self):
        super(Akeypre, self).__init__()
        self.allowed_values = ["unif", "geom"]
    @property
    def name(self):
        return "akeypre"

class Akeysolve(OneStringParameter):
    """ Usage: akeySOLVE {key}

    Specifies how the the finite element mesh should be adaptively subdivided during the
    solve-estimate-refine iterations of a fe-manual finite element calculation. This allows
    for various a posteriori refinement schemes.  key is a text string that specifies the method
    used to guide adaptive refinement:

    * resi - Residual-based a posteriori refinement
    """
    def __init__(self):
        super(Akeysolve, self).__init__()
        self.allowed_values = ["resi"]
    @property
    def name(self):
        return "akeysolve"

class Async(OneIntegerParameter):
    """ This optional keyword allows users to perform the different tasks in a mg-para parallel
    run asynchronously. Specifically, a processor masquerades as process rank in a parallel
    focusing run and provides output (data files and energies/forces) appropriate to that
    processor's local partition. The user must then assemble the results after all processes
    complete. First, this option is useful for scheduling on-demand resources: this makes it
    easy for users to backfill into the available processes in a queue. Second, this option is
    useful for running on limited resources: this enables users without access to large parallel
    machines to still perform the same calculations. The syntax is

    async { rank }

    where rank is the integer ID of the particular processor to masquerade as. Processor IDs
    range from 0 to N-1, where N is the total number of processors in the run (see pdime).
    Processor IDs are related to their position in the overall grid by p = nx ny k + nx j + i
    where nx is the number of processors in the x-direction, ny is the number of processors in
    the y-direction, nz is the number of processors in the z-direction, i is the index of the
    processor in the x-direction, j is the index of the processor in the y-direction, k is the
    index of the processor in the z-direction, and p is the overall rank of the processor."""
    @property
    def name(self):
        return "async"

class Bcfl(OneStringParameter):
    """ Specifies the type of boundary conditions used to solve the Poisson-Boltzmann equation.
    The syntax is

    bcfl {flag}

    where flag is a text string that identifies the type of conditions to be used:

    * zero - "Zero" boundary condition. Dirichlet conditions where the potential at the boundary
    is set to zero. This condition is not commonly used and can result in large errors if used
    inappropriately.
    * sdh - "Single Debye-Huckel" boundary condition. Dirichlet condition where the potential at
    the boundary is set to the values prescribed by a Debye-Huckel model for a single sphere with
    a point charge, dipole, and quadrupole. The sphere radius in this model is set to the radius
    of the biomolecule and the sphere charge, dipole, and quadrupole are set to the total moments
    of the protein. This condition works best when the boundary is sufficiently far from the
    biomolecule.
    * mdh - "Multiple Debye-Huckel" boundary condition. Dirichlet condition where the potential
    at the boundary is set to the values prescribed by a Debye-Huckel model for a multiple,
    non-interacting spheres with a point charges. The radii of the non-interacting spheres are
    set to the atomic radii of and the sphere charges are set to the atomic charges. This
    condition works better than sdh for closer boundaries but can be very slow for large
    biomolecules.
    * focus - "Focusing" boundary condition. Dirichlet condition where the potential at the
    boundary is set to the values computed by the previous (usually lower-resolution) PB
    calculation. This is used in sequential focusing performed manually in mg-manual
    calculations. All of the boundary points should lie within the domain of the previous
    calculation for best accuracy; if any boundary points lie outside, their values are computed
    using single Debye-Huckel boundary conditions (see above).
    * map - Specifying map allows a previously calculated potential map to be used in a new
    focusing calculation. A typical scenario is using the same coarse grid for multiple focusing
    calculations. A potential map can be written once from a coarse grid calculation, then used
    in subsequent runs to bypass the need to recalculate the coarse grid. See the READ keyword
    pot and the attached example files for its use.  NOTE:  this functionality is only available
    in the current developmental release of APBS."""
    def __init__(self):
        super(Bcfl, self).__init__()
        self.allowed_values = ["zero", "sdh", "mdh", "focus", "map"]
    @property
    def name(self):
        return "bcfl"

class Gcent(Parameter):
    """ This class serves as both the storage class for the gcent keyword _and_ a class from which
    cgcent and fgcent are inherited.

    Specify the center of the grid based on a molecule's center or absolute coordinates for a
    mg-manual multigrid calculation. The syntax is:

    gcent { mol id | xcent ycent zcent }

    where the user can specify either:

    * mol {id} - Center the grid on molecule with integer ID id; as assigned in the READ section.
    Molecule IDs are assigned in the order they are read, starting at 1.

    or the user can specify

    * xcent ycent zcent - The floating point coordinates (in &Aring;) at which the grid is centered.
    Based on the PDB coordinate frame."""
    def __init__(self):
        self.content_dict = {"mol" : None, "xcent" : None, "ycent" : None, "zcent" : None}
    @property
    def name(self):
        return "gcent"
    def parse(self, tokens):
        token = tokens.pop(0).lower()
        if token == "mol":
            id_token = int(tokens.pop(0))
            self.content_dict["mol"] = id_token
        else:
            self.content_dict["xcent"] = float(token)
            self.content_dict["ycent"] = float(tokens.pop(0))
            self.content_dict["zcent"] = float(tokens.pop(0))
    def validate(self):
        # No validation here since it happens in the parsing above
        pass
    def contents(self):
        return self.content_dict
    def __str__(self):
        outstr = "%s " % self.name
        mol_id = self.content_dict["mol"]
        if mol_id:
            outstr = outstr + "mol %d" % mol_id
        else:
            outstr = outstr + "%g %g %g" % (self.content_dict["xcent"], self.content_dict["ycent"],
                                            self.content_dict["zcent"])
        return outstr

class Cgcent(Gcent):
    """ Specify the center of the coarse grid (in a focusing calculation) based on a molecule's
    center or absolute coordinates for a multigrid (mg-manual, mg-auto, mg-para) Poisson-Boltzmann
    calculation. The syntax is

    cgcent { mol id | xcent ycent zcent }

    The arguments for this keyword are either
    * mol id - Center the grid on molecule with integer ID id; as assigned in the READ section with
    a READ mol command.

    or

    * xcent ycent zcent - Center the grid on the (floating point) coordinates (in &Aring;) at which
    the grid is centered. Based on the PDB coordinate frame.

    See also: gcent, fgcent """
    @property
    def name(self):
        return "cgcent"

class Cglen(Glen):
    """ Specify the length of the coarse grid (in a focusing calculation) for an automatic multigrid
    (mg-auto, mg-para) Poisson-Boltzmann calculation.  This may be different in each direction. Its
    syntax is

    cglen {xlen ylen zlen}

    This is the starting mesh, so it should be large enough to complete enclose the biomolecule and
    ensure that the chosen boundary condition (see bcfl) is appropriate.

    xlen ylen zlen - Grid lengths (floating point numbers) in the x-, y-, and z-directions in
    &Aring.

    See also: fglen or glen """
    @property
    def name(self):
        return "cglen"

class Chgm(OneStringParameter):
    """ Specify the method by which the biomolecular point charges (i.e., Dirac delta functions) by
    which charges are mapped to the grid for a multigrid (mg-manual, mg-auto, mg-para) Poisson-
    Boltzmann calculation.  As we are attempting to model delta functions, the support (domain) of
    these discretized charge distributions is always a function of the grid spacing. The syntax for
    this command is:

    chgm {flag}

    where flag is a text string that specifies the type of discretization:

    * spl0 - Traditional trilinear interpolation (linear splines). The charge is mapped onto the
    nearest-neighbor grid points. Resulting potentials are very sensitive to grid spacing, length,
    and position.
    * spl2 - Cubic B-spline discretization. The charge is mapped onto the nearest- and next-nearest-
    neighbor grid points. Resulting potentials are somewhat less sensitive (than spl0) to grid
    spacing, length, and position.
    * spl4 - Quintic B-spline discretization. Similar to spl2, except the charge/multipole is
    additionally mapped to include next-next-nearest neighbors (125 grid points receive charge
    density). """
    def __init__(self):
        super(Chgm, self).__init__()
        self.allowed_values = ["spl0", "spl2", "spl4"]
    @property
    def name(self):
        return "chgm"

class Dime(ThreeIntegerParameter):
    """ Specifies the number of grid points per processor for grid-based discretization. Its syntax
    is

    dime {nx ny nz}

    For mg-manual calculations, the arguments are dependent on the choice of nlev by the formula

    n = c 2^{l + 1} + 1

    where n is the dime argument, c is a non-zero integer, l is the nlev value. The most common
    values for grid dimensions are 65, 97, 129, and 161 (they can be different in each direction);
    these are all compatible with a nlev value of 4. If you happen to pick a "bad" value for the
    dimensions (i.e., mismatch with nlev), the APBS code will adjust the specified dime downwards
    to more appropriate values. This means that "bad" values will typically result in lower
    resolution/accuracy calculations! The arguments for this keyword are:

    nx ny nz

    the (integer) number of grid points in the x-, y-, and z-directions, respectively.

    NOTE: dime should be interpreted as the number of grid points per processor for all
    calculations, including mg-para. This interpretation helps manage the amount of memory
    per-processor -- generally the limiting resource for most calculations."""
    @property
    def name(self):
        return "dime"
    def fix_dimension(self, in_dim):
        """ Check the dimension for compatibility with multigrid.  Return the same dimension or a
        fixed value if there are problems. """
        facs = factors(in_dim-1)
        nlev = facs.count(2)-1
        while nlev < 2:
            facs.sort()
            facs.pop(0)
            facs.append(2)
            nlev = nlev + 1
        out_dim = product(facs)
        return out_dim+1
    def validate(self):
        newval = self.fix_dimension(self.xint)
        if newval != self.xint:
            _LOGGER.error("Corrected dimension %d to %d for compatibility with \
multigrid.\n", self.xint, newval)
            self.xint = newval
        newval = self.fix_dimension(self.yint)
        if newval != self.yint:
            _LOGGER.error("Corrected dimension %d to %d for compatibility with \
multigrid.\n", self.yint, newval)
            self.yint = newval
        newval = self.fix_dimension(self.zint)
        if newval != self.zint:
            _LOGGER.error("Corrected dimension %d to %d for compatibility with \
multigrid.\n", self.zint, newval)
            self.zint = newval

class Domainlength(ThreeFloatParameter):
    """ Specify the rectangular finite element mesh domain lengths for fe-manual finite element
    calculations.  This length may be different in each direction. If the usemesh keyword is
    included, then this command is ignored. The syntax is:

    domainLength {xlen ylen zlen}

    where the parameters xlen, ylen, zlen are floating point numbers that specify the mesh lengths
    in the x-, y-, and z-directions (respectively) in units of &Aring;. """
    @property
    def name(self):
        return "domainlength"
    def validate(self):
        if (self.xfloat < FLOAT_EPSILON) or (self.yfloat < FLOAT_EPSILON) or (self.zfloat
                                                                              < FLOAT_EPSILON):
            errstr = "One of the domain length parameters is zero or negative"
            raise ValueError(errstr)

class Ekey(OneStringParameter):
    """ Specify the method used to determine the error tolerance in the solve-estimate-refine
    iterations of the finite element solver (fe-manual). The syntax is:

    ekey { flag }

    where flag is a text string that determines the method for error calculation:

    * simp - Per-simplex error limit
    * global - Global (whole domain) error limit
    * frac - Fraction of simplices you'd like to see refined at each iteration

    See also: etol """
    def __init__(self):
        super(Ekey, self).__init__()
        self.allowed_values = ["simp", "global", "frac"]
    @property
    def name(self):
        return "ekey"

class Etol(OneFloatParameter):
    """ Error tolerance

    Finite difference multigrid methods

    Current developmental releases of APBS provide support for the optional etol keyword to specify
    the tolerance for iterations of the PMG partial differential equation solver.  The syntax is

    etol { tol }

    where tol is the (floating point) numerical value for the error tolerance.  This keyword is
    optional and is intended for mg-manual, mg-auto, and mg-para calculation types.

    Finite element methods

    Specify the tolerance for error-based adaptive refinement during the solve-estimate-refine
    iterations of the finite element solver (fe-manual). The syntax is

    etol { tol }

    where tol is the (floating point) numerical value for the error tolerance.

    See also: ekey """
    @property
    def name(self):
        return "etol"

class Fgcent(Gcent):
    """ Specify the center of the fine grid (in a focusing calculation) based on a molecule's center
    or absolute coordinates for mg-para and mg-auto multigrid calculations. The syntax is:

    fgcent { mol id | xcent ycent zcent }

    where a user can specify either

    * mol {id} - Center the grid on molecule with integer ID id; as assigned in the READ section of
    the input file. Molecule IDs are assigned in the order they are read, starting at 1.

    or the user can specify

    * xcent ycent zcent - Center the grids on the coordinates (floating point numbers in &Aring;) at
    which the grid is centered. Based on the input molecule PDB coordinate frame.

    See also: cgcent """
    @property
    def name(self):
        return "fgcent"

class Fglen(Glen):
    """ Specifies the fine mesh domain lengths in a multigrid focusing calculation (mg-para or
    mg-auto); this may be different in each direction. The syntax is

    fglen {xlen ylen zlen}

    This should enclose the region of interest in the molecule. The arguments to this command are:

    xlen ylen zlen - Grid lengths (floating point numbers) in the x-, y-, and z-directions in
    &Aring;.

    See also: cglen """
    @property
    def name(self):
        return "fglen"

class Grid(ThreeFloatParameter):
    """ Specify the mesh grid spacings for multigrid mg-manual calculations.  This value may be
    different in each direction. The syntax is:

    grid {hx hy hz}

    where hx hy hz are the (floating point) grid spacings in the x-, y-, and z-directions
    (respectively) in &Aring;. """
    @property
    def name(self):
        return "grid"
    def validate(self):
        if (self.xfloat < FLOAT_EPSILON) or (self.yfloat < FLOAT_EPSILON) or (self.zfloat
                                                                              < FLOAT_EPSILON):
            errstr = "All grid spacings must be greater than zero!"
            raise ValueError(errstr)

class Ion(Parameter):
    """ Specify the bulk concentrations of mobile ion species present in the system. This command
    can be repeated as necessary to specify multiple types of ions; however, only the largest ionic
    radius is used to determine the ion-accessibility function. The total bulk system of ions must
    be electroneutral which means the charge densities/concentrations of positive and negative ions
    must be equal. The syntax is

    ion charge {charge} conc {conc} radius {radius}

    where

    * charge - Mobile ion species charge (floating point number in ec)
    * conc - Mobile ion species concentration (floating point number in M)
    * radius - Mobile ion species radius (floating point number in &Aring;) """
    def __init__(self):
        self.charge = None
        self.conc = None
        self.radius = None
    @property
    def name(self):
        return "ion"
    def validate(self):
        if self.conc < 0:
            errstr = "Concentration is negative (%g)" % self.conc
            raise ValueError(errstr)
        if self.radius < FLOAT_EPSILON:
            errstr = "Radius is zero or negative (%g)" % self.radius
            raise ValueError(errstr)
    def raise_error(self, token):
        """ General method for raising errors about bad keyword values """
        errstr = "Unexpected token (%s) in ION keyword" % token
        raise ValueError(errstr)
    def parse(self, tokens):
        token = tokens.pop(0).lower()
        if token != "charge":
            self.raise_error(token)
        self.charge = int(tokens.pop(0))
        token = tokens.pop(0).lower()
        if token != "conc":
            self.raise_error(token)
        self.conc = float(tokens.pop(0))
        token = tokens.pop(0).lower()
        if token != "radius":
            self.raise_error(token)
        self.radius = float(tokens.pop(0))
    def contents(self):
        return {"charge" : self.charge, "conc" : self.conc, "radius" : self.radius}
    def __str__(self):
        outstr = "ion charge %d conc %g radius %g" % (self.charge, self.conc, self.radius)
        return outstr

class Eqntype(OneStringParameter):
    """ Type of equation being solved.

    eqntype {type}

    where type is one of

    * lpbe - Linearized Poisson-Boltzmann
    * lrpbe - Linearized regularized Poisson-Boltzmann
    * npbe - Nonlinear Poisson-Boltzmann
    * nrpbe - Nonlinear regularized Poisson-Boltzmann """
    # TODO - This isn't documented yet; needs to be updated in documentation, etc.
    def __init__(self):
        super(Eqntype, self).__init__()
        self.allowed_values = ["lpbe", "lrpbe", "npbe", "nrpbe"]
    @property
    def name(self):
        return "eqntype"

class Maxsolve(OneIntegerParameter):
    """ Specify the number of times to perform the solve-estimate-refine iteration of the finite
    element solver (fe-manual). The syntax is

    maxsolve { num }

    where num is an integer indicating the desired maximum number of iterations.

    See also: maxvert, targetRes  """
    @property
    def name(self):
        return "maxsolve"
    def validate(self):
        if self.parm < 1:
            errstr = "maxsolve is less than 1"
            raise ValueError(errstr)

class Maxvert(OneIntegerParameter):
    """ Specify the maximum number of vertices to allow during solve-estimate-refine cycle of
    finite element solver (fe-manual). This places a limit on the memory that can be used by the
    solver. The syntax is

    maxvert { num }

    where num is an integer indicating the maximum number of vertices.

    See also: targetRes, maxsolve. """
    @property
    def name(self):
        return "maxvert"
    def validate(self):
        if self.parm < 1:
            raise ValueError("maxvert is less than 1")

class Mol(OneIntegerParameter):
    """ Specify the molecule for which the PBE is to be solved. IDs are based on the order in which
    molecules are read by READ mol statements, starting from 1. The syntax is

    mol {id}

    where id is the integer ID of the molecule for which the Poisson-Boltzmann equation is to be
    solved. """
    @property
    def name(self):
        return "mol"

class Nlev(OneIntegerParameter):
    """ Specify the depth of the multilevel hierarchy used in the mg-manual multigrid solver. See
    dime for a discussion of how nlev relates to grid dimensions. The syntax is

    nlev {lev}

    where lev is an integer indicating the desired depth of the multigrid hierarchy. """
    @property
    def name(self):
        return "nlev"
    def validate(self):
        if self.parm < 1:
            raise ValueError("nlev is less than 1")

class Ofrac(OneFloatParameter):
    """ Specify the amount of overlap to include between the individual processors meshes in a
    parallel focusing calculation (mg-para). This should be a value between 0 and 1. The syntax is

    ofrac {frac}

    where frac is a floating point value between 0.0 and 1.0 denoting the amount of overlap between
    processors.  Empirical evidence suggests that an ofrac value of 0.1 is sufficient to generate
    stable energies. However, this value may not be sufficient to generate stable forces and/or good
    quality isocontours. """
    @property
    def name(self):
        return "ofrac"
    def validate(self):
        if self.parm < FLOAT_EPSILON:
            raise ValueError("ofrac is less than or equal to 0")

class Pdie(OneFloatParameter):
    """ Specify the dielectric constant of the biomolecule. This is usually a value between 2 to 20,
    where lower values consider only electronic polarization and higher values consider additional
    polarization due to intramolecular motion. The syntax is:

    pdie {diel}

    where diel is the floating point value of the unitless biomolecular dielectric constant.

    See also: sdie """
    @property
    def name(self):
        return "pdie"
    def validate(self):
        if self.parm < FLOAT_EPSILON:
            raise ValueError("pdie is zero or negative")

class Pdime(ThreeIntegerParameter):
    """ Specify the processor array to be used in a parallel focusing (mg-para) calculation. The
    product npx * npy * npz should be less than or equal to the total number of processors with
    which APBS was invoked (usually via mpirun). If more processors are provided at invocation than
    actually used during the run, the extra processors are not used in the calculation. The
    processors are tiled across the domain in a Cartesian fashion with a specified amount of overlap
    (see ofrac) between each processor to ensure continuity of the solution. Each processor's
    subdomain will contain the number of grid points specified by the dime keyword. The syntax is:

    pdime {npx npy npz}

    where npx npy npz are the integer number of processors to be used in the x-, y- and z-directions
    of the system. For broad spatial support of the splines, every charge included in partition
    needs to be at least 1 grid space (chgm spl0), 2 grid spaces (chgm spl2), or 3 grid spaces (chgm
    spl4) away from the partition boundary.

    See also: ofrac """
    @property
    def name(self):
        return "pdime"

    def validate(self):
        if (self.xint < 1) or (self.yint < 1) or (self.zint < 1):
            errstr = "One of the dimensions is less than 1"
            raise ValueError(errstr)

class Sdens(OneFloatParameter):
    """ Specify the number of grid points per square-angstrom to use in discontinuous surface
    constructions (e.g., molecular surface and solvent-accessible surfaces). Ignored when srad is
    0.0 or srfm is spl2. There is a direct correlation between this value used for the surface
    sphere density, the accuracy of the surface calculations, and the APBS calculation time. The
    APBS "suggested" value is 10.0. The syntax of this command is

    sdens {density}

    where density is the floating point surface sphere density (in grid points/&Aring;^2).

    See also: srfm """
    @property
    def name(self):
        return "sdens"
    def validate(self):
        if self.parm < FLOAT_EPSILON:
            raise ValueError("sdens is zero or negative")

class Sdie(OneFloatParameter):
    """ Specify the dielectric constant of the solvent. Bulk water at biologically-relevant
    temperatures is usually modeled with a dielectric constant of 78-80. The syntax is

    sdie {diel}

    where diel is a floating point number representing the solvent dielectric constant (unitless).

    See also: pdie """
    @property
    def name(self):
        return "sdie"
    def validate(self):
        if self.parm < FLOAT_EPSILON:
            raise ValueError("sdie is zero or negative")

class Srad(OneFloatParameter):
    """ Specify the radius of the solvent molecules; this parameter is used to define the dielectric
    function for probe-based dielectric definitions (see srfm). This value is usually set to 1.4
    &Aring; for water. This keyword is ignored when any of the spline-based surfaces are used (e.g.,
    spl2, see srfm), since they are not probe-based. The syntax for this command is:

    srad {radius}

    where radius is the floating point solvent radius (in &Aring).

    See also: srfm """
    @property
    def name(self):
        return "srad"
    def validate(self):
        if self.parm < 0:
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
    def __init__(self):
        super(Srfm, self).__init__()
        self.allowed_values = ["mol", "smol", "spl2", "spl4", "sacc"]
    @property
    def name(self):
        return "srfm"

class Swin(OneFloatParameter):
    """ Specify the size of the support (i.e., the rate of change) for spline-based surface
    definitions (see srfm). Usually 0.3 &Aring;. The syntax is:

    swin {win}

    where win is a floating point number for the spline window width (in &Aring;). Note that, per
    the analysis of Nina, Im, and Roux (doi:10.1016/S0301-4622(98)00236-1), the force field
    parameters (radii) generally need to be adjusted if the spline window is changed. """
    @property
    def name(self):
        return "swin"
    def validate(self):
        if self.parm < FLOAT_EPSILON:
            raise ValueError("swin is zero or negative")

class Targetnum(OneIntegerParameter):
    """ Specify the target number of vertices in the initial finite element mesh for fe-manual
    calculations.  Initial refinement will continue until this number is reached or the the longest
    edge of every simplex is below targetNum. The syntax is

    targetNum { num }

    where num is an integer denoting the target number of vertices in initial mesh. See also:
    targetRes """
    @property
    def name(self):
        return "targetnum"
    def validate(self):
        if self.parm < 1:
            raise ValueError("targetnum is less than 1")

class Targetres(OneFloatParameter):
    """ Specify the target resolution of the simplices in a finite element mesh (fe-manual);
    refinement will continue until the longest edge of every simplex is below this value. The syntax
    is

    targetRes { res }

    where res is a floating point number denoting the target resolution for longest edges of
    simplices in mesh (in &Aring;).

    See also: maxvert, maxsolve, targetNum. """
    @property
    def name(self):
        return "targetres"
    def validate(self):
        if self.parm < FLOAT_EPSILON:
            raise ValueError("targetres is zero or negative")

class Usemap(Parameter):
    """ Specify pre-calculated coefficient maps to be used in the Poisson-Boltzmann calculation.
    These must have been input via an earlier READ statement. The syntax is

    usemap {type} {id}

    where the mandatory keywords are

    * type - A string that specifies the type of pre-calculated map to be read in:
      - diel - Dielectric function map (as read by read diel); this causes the pdie, sdie, srad,
      swin, and srfm parameters and the radii of the biomolecular atoms to be ignored when computing
      dielectric maps for the Poisson-Boltzmann equation. Note that the pdie and sdie values are
      still used for some boundary condition calculations as specified by bcfl.
      - kappa - Mobile ion-accessibility function map (as read by read kappa); this causes the swin
      and srfm parameters and the radii of the biomolecular atoms to be ignored when computing
      mobile ion values for the Poisson-Boltzmann equation. The ion parameter is not ignored and
      will still be used.
      - charge - Charge distribution map (as read by read charge); this causes the chgm parameter
      and the charges of the biomolecular atoms to be ignored when assembling the fixed charge
      distribution for the Poisson-Boltzmann equation.
      - pot - Potential map (as read by read pot); this option requires setting bcfl to map.  NOTE:
      this functionality is only available in the current developmental release of APBS.
    * id - As described in the READ command documentation, this integer ID specifies the particular
    map read in with READ. These IDs are assigned sequentially, starting from 1, and incremented
    independently for each map type read by APBS. In other words, a calculation that uses two PQR
    files, one parameter file, three charge maps, and four dielectric maps would have PQR files
    with IDs 1-2, a parameter file with ID 1, charge maps with IDs 1-3, and dielectric maps with
    IDs 1-4. """
    def __init__(self):
        super(Usemap, self).__init__()
        self.type = None
        self.allowed_values = ["diel", "kappa", "charge", "pot"]
        self.map_id = None
    @property
    def name(self):
        return "usemap"
    def parse(self, tokens):
        token = tokens.pop(0).lower()
        if token in self.allowed_values:
            self.type = token
        else:
            errstr = "Unknown map type (%s) in usemap" % token
            raise ValueError(errstr)
        self.map_id = int(tokens.pop(0))
    def validate(self):
        # Validation happens in parsing
        pass
    def contents(self):
        return {"type" : self.type, "map_id" : self.map_id}
    def __str__(self):
        outstr = "usemap %s %d" % (self.type, self.map_id)
        return outstr

class Usemesh(OneIntegerParameter):
    """ Specify the external finite element mesh to be used in the finite element Poisson-Boltzmann
    calculation (fe-manual). These must have been input via an earlier READ mesh statement. The
    syntax is

    usemesh {id}

    where id is an integer ID specifying the particular map read in with READ mesh. These IDs are
    assigned sequentially, starting from 1, and incremented independently for each mesh read by
    APBS. """
    @property
    def name(self):
        return "usemesh"

class Write(Parameter):
    """ This controls the output of scalar data calculated during the Poisson-Boltzmann run. This
    keyword can be repeated several times to provide various types of data output from APBS. The
    syntax is

    write {type} {format} {stem}

    * type - A string indicating what type of data to output:
      - charge - Write out the biomolecular charge distribution in units of ec (electron charge)
      per &Aring;^3. (multigrid only)
      - pot - Write out the electrostatic potential in units of kb T ec-1. (multigrid and finite
      element)
      - atompot - Write out the electrostatic potential at each atom in units of kb T ec^-1. The
      format for the output must be specified as flat. The values are listed sequentially from 1 to
      NATOMS. See below for more information about format types. (multigrid and finite element).
      - smol - Write out the solvent accessibility defined by the molecular surface definition (see
      srfm smol). Values are unitless and range from 0 (inaccessible) to 1 (accessible). (multigrid
      and finite element)
      - sspl - Write out the spline-based solvent accessibility (see srfm spl2). Values are unitless
      and range from 0 (inaccessible) to 1 (accessible) (multigrid and finite element)
      - vdw - Write out the van der Waals-based solvent accessibility (see srfm smol with srad 0.0).
      Values are unitless and range from 0 (inaccessible) to 1 (accessible). (multigrid and finite
      element)
      - ivdw - Write out the inflated van der Waals-based ion accessibility (see srfm smol). Values
      are unitless and range from 0 (inaccessible) to 1 (accessible). (multigrid and finite element)
      - lap - Write out the Laplacian of the potential in units of kB T ec^-1 &Aring;^-2. (multigrid
      only)
      - edens - Write out the "energy density" in units of kB T ec^-1 &Aring;^-2. (multigrid only)
      - ndens - Write out the total mobile ion number density for all ion species in units of M.
      (multigrid only)
      - qdens - Write out the total mobile charge density for all ion species in units of ec M.
      (multigrid only)
      - dielx - Write out the dielectric map shifted by 1/2 grid spacing in the x-direction (see
      READ diel). The values are unitless. (multigrid only)
      - diely - Write out the dielectric map shifted by 1/2 grid spacing in the y-direction (see
      READ diel). The values are unitless. (multigrid only)
      - dielz - Write out the dielectric map shifted by 1/2 grid spacing in the z-direction (see
      READ diel). The values are unitless. (multigrid only)
      - kappa - Write out the ion-accessibility kappa map (see READ kappa). The values are in units
      of &Aring;^-2 (multigrid only)
    * format - A string that specifies the format for writing out the data.
      - dx - Write out data in OpenDX format. This is the preferred format for APBS I/O. (multigrid
      and finite element).
      - avs - Write out data in AVS UCD format. (finite element only)
      - uhbd - Write out data in UHBD format. (multigrid only)
      - gz - Write out OpenDX data in gzipped (zlib) compatible format. Appends .dx.gz to the
      filename.
      - flat - Write out data as a plain text file. (multigrid and finite element).
    * stem - A string that specifies the path for the output; files are written to stem.XYZ, where
    XYZ is determined by the file format (and processor rank for parallel calculations). If the
    pathname contains spaces, then it must be surrounded by double quotes; e.g., \"/path with
    spaces/foo.in\". """
    def __init__(self):
        super(Write, self).__init__()
        self.allowed_type_values = ["charge", "pot", "atompot", "smol", "sspl", "vdw", "ivdw",
                                    "lap", "edens", "ndens", "qdens", "dielx", "diely", "dielz",
                                    "kappa"]
        self.allowed_format_values = ["avs", "uhbd", "gz", "flat", "dx"]
        self.type = None
        self.format = None
        self.stem = None
    @property
    def name(self):
        return "write"
    def parse(self, tokens):
        self.type = tokens.pop(0).lower()
        self.format = tokens.pop(0).lower()
        self.stem = tokens.pop(0)
    def validate(self):
        if not self.type in self.allowed_type_values:
            errstr = "Unexpected type token (%s) for write" % self.type
            raise ValueError(errstr)
        if not self.format in self.allowed_format_values:
            errstr = "Unexpected format token (%s) for write" % self.format
            raise ValueError(errstr)
    def contents(self):
        return {"type" : self.type, "format" : self.format, "stem" : self.stem}
    def __str__(self):
        return "write %s %s \"%s\"" % (self.type, self.format, self.stem)

class Writemat(FormatPathParameter):
    """ This controls the output of the mathematical operators in the Poisson-Boltzmann equation as
    matrices in Harwell-Boeing matrix format (multigrid only). The syntax is:

    writemat {type} {stem}

    where

    * type - A string that indicates what type of operator to output:
      - poisson - Write out the Poisson operator - nabla cdot epsilon nabla.
    *stem - A string that specifies the path for the matrix """
    def __init__(self):
        super(Writemat, self).__init__()
        self.type = None
        self.stem = None
        self.allowed_values = ["poisson"]
    @property
    def name(self):
        return "writemat"
    def parse(self, tokens):
        self.type = tokens.pop(0).lower()
        self.stem = tokens.pop(0)
    def validate(self):
        if not self.type in self.allowed_values:
            errstr = "Unknown token (%s) for %s" % self.type
            raise ValueError(errstr)
    def contents(self):
        return {"type" : self.type, "stem" : self.stem}
    def __str__(self):
        return "writemat %s %s\n" % (self.type, self.stem)

class Elec(ParameterSection):
    """ The ELEC block of an APBS input file is used for polar solvation (electrostatics)
    calculations and has the following syntax:

    ELEC [ name {id} ]
        {type}
        {keywords...}
    END

    where the indentation and linefeeds are included for clarity; only whitespace is needed in the
    input file.  The {id} tag allows the user to name ELEC blocks, as described in the ELEC block
    naming section.  The {type} command defines the Types of ELEC calculation to be performed.
    Finally, the {keywords} are calculation-specific commands that customize the particular type of
    calculation.

    This section is the main component for polar solvation calculations in APBS runs. There may be
    several ELEC sections, operating on different molecules or using different parameters for
    multiple runs on the same molecule. The order of the ELEC statement can matter since certain
    types of boundary conditions (bcfl) can require information about previous calculations. """
    def __init__(self):
        super(Elec, self).__init__()
        self.my_name = "elec"
    def parse(self, tokens):
        token = tokens.pop(0)
        while True:
            token_name = token.lower()
            if token_name in ["name", "solvtype", "akeypre", "akeysolve", "async", "bcfl",
                              "calcenergy", "calcforce", "cgcent", "cglen", "chgm", "dime",
                              "domainlength", "ekey", "etol", "fgcent", "fglen", "gcent", "glen",
                              "grid", "eqntype", "maxsolve", "maxvert", "mol", "nlev", "ofrac",
                              "pdie", "pdime", "sdens", "sdie", "srad", "srfm", "swin", "targetnum",
                              "targetres", "temp", "usemap", "usemesh", "write", "writemat"]:
                # This section handles simple command statements
                self.create_store_single_object(token_name, tokens)
            elif token_name in ["ion"]:
                # This section handles command statements that can appear multiple times
                self.create_store_multiple_objects(token_name, tokens)
            elif token_name in Solvtype().allowed_values:
                # This is a special section to handle the old-format ELEC solver type declaration
                solvtype = Solvtype()
                solvtype.parse([token_name])
                solvtype.validate()
                self.content_dict["solvtype"] = solvtype
            elif token_name in Eqntype().allowed_values:
                # This is a special section to handle the old-format ELEC equation type declaration
                eqntype = Eqntype()
                eqntype.parse([token_name])
                eqntype.validate()
                self.content_dict["eqntype"] = eqntype
            elif token_name == "end":
                return
            else:
                errstr = "Unknown token (%s) for %s" % (token, self.name)
                raise ValueError(errstr)
            token = tokens.pop(0)
    def validate_mgauto(self):
        """ Validate current input file contents """
        required_parms = ["bcfl", "chgm", "cgcent", "cglen", "dime", "fgcent", "fglen",
                          "eqntype", "mol", "pdie", "sdens", "sdie", "srad", "srfm", "swin",
                          "temp"]
        for parm in required_parms:
            if not parm in self.content_dict:
                errstr = "Missing required parameter %s for %s" % (parm,
                                                                   self.content_dict["solvtype"])
                raise ValueError(errstr)
    def validate_mgpara(self):
        """ Validate current input file contents """
        required_parms = ["bcfl", "chgm", "cgcent", "cglen", "dime", "fgcent", "fglen",
                          "eqntype", "mol", "ofrac", "pdie", "pdime", "sdens", "sdie", "srad",
                          "srfm", "swin", "temp"]
        for parm in required_parms:
            if not parm in self.content_dict:
                errstr = "Missing required parameter %s for %s" % (parm,
                                                                   self.content_dict["solvtype"])
                raise ValueError(errstr)
    def validate_mgmanual(self):
        """ Validate current input file contents """
        required_parms = ["bcfl", "chgm", "dime", "gcent", "eqntype", "mol", "nlev", "pdie",
                          "sdens", "sdie", "srad", "srfm", "swin", "temp"]
        for parm in required_parms:
            if not parm in self.content_dict:
                errstr = "Missing required parameter %s for %s" % (parm,
                                                                   self.content_dict["solvtype"])
                raise ValueError(errstr)
        if (not "glen" in self.content_dict) and (not "grid" in self.content_dict):
            errstr = "Either grid or glen needs to be specified for mg-manual"
            raise ValueError(errstr)
        if self.content_dict["solvtype"] == "mg-manual":
            dime = self.content_dict["dime"]
            nlev = self.content_dict["nlev"]
            ncurr = nlev.parm
            for dim in [dime.xint, dime.yint, dime.zint]:
                facs = factors(dim-1)
                ncalc = facs.count(2)-1
                if ncalc < ncurr:
                    errstr = "Your given nlev (%d) is larger than the permissible value (%d)" \
                             % (ncurr, ncalc)
                    raise ValueError(errstr)
    def validate_femanual(self):
        """ Validate current input file contents """
        required_parms = ["akeypre", "akeysolve", "bcfl", "chgm", "domainlength", "ekey",
                          "etol", "eqntype", "maxsolve", "maxvert", "mol", "pdie", "sdens",
                          "sdie", "srad", "srfm", "swin", "temp"]
        for parm in required_parms:
            if not parm in self.content_dict:
                errstr = "Missing required parameter %s for %s" % (parm,
                                                                   self.content_dict["solvtype"])
                raise ValueError(errstr)
    def validate(self):
        """ Validate current input file contents """
        solvtype = self.content_dict["solvtype"].parm
        if solvtype == "mg-auto":
            self.validate_mgauto()
        elif solvtype == "mg-para":
            self.validate_mgpara()
        elif solvtype in ["mg-manual", "mg-dummy"]:
            self.validate_mgmanual()
        elif solvtype == "fe-manual":
            self.validate_femanual()
        else:
            raise TypeError("Don't know how to process solver type %s" % solvtype)
    def __str__(self):
        outstr = self.format_block_start()
        keys = list(self.content_dict.keys())
        keys.remove("name")
        outstr = outstr + "\t%s\n" % self.content_dict["solvtype"]
        keys.remove("solvtype")
        outstr = outstr + "\t%s\n" % self.content_dict["eqntype"]
        keys.remove("eqntype")
        for key in keys:
            values = self.content_dict[key]
            try:
                for value in values:
                    outstr = outstr + "\t%s\n" % value
            except TypeError:
                outstr = outstr + "\t%s\n" % values
        outstr = outstr + "end\n"
        return outstr

### These are the APOLAR-specific classes

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

class Apolar(ParameterSection):
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
        self.my_name = "apolar"
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
            if parm not in self.content_dict:
                errstr = "APOLAR missing required parameter %s" % parm
                _LOGGER.error("Incomplete parameters: %s", self.content_dict)
                raise ValueError(errstr)
    def __str__(self):
        outstr = self.format_block_start()
        keys = list(self.content_dict.keys())
        keys.remove("name")
        for key in keys:
            values = self.content_dict[key]
            try:
                for value in values:
                    outstr = outstr + "\t%s\n" % value
            except TypeError:
                outstr = outstr + "\t%s\n" % values
        outstr = outstr + "end\n"
        return outstr
