""" Handle the storage of APBS ELEC block input file parameters """
import sys
import logging

_LOGGER = logging.getLogger("elec-parser")

from .utility import factors, product
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

class Glen(parameter.ThreeFloatParameter):
    """ Specify the mesh domain lengths for multigrid mg-manual calculations.  These lengths may be
    different in each direction. The syntax is:

    glen {xlen ylen zlen}

    where xlen ylen zlen are the (floating point) grid lengths in the x-, y-, and z-directions
    (respectively) in &Aring;.

    See also: grid  """
    def __init__(self):
        super(Glen, self).__init__()
        self._short_name = "glen"
    def validate(self):
        if (self.xfloat < parameter.FLOAT_EPSILON) \
        or (self.yfloat < parameter.FLOAT_EPSILON) or (self.zfloat < parameter.FLOAT_EPSILON):
            errstr = "One of the grid lengths is zero or negative (%g, %g, %g)" % \
            (self.xfloat, self.yfloat, self.zfloat)
            raise ValueError(errstr)

class Solvtype(parameter.OneStringParameter):
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
        self._allowed_values = ["mg-auto", "mg-para", "mg-manual", "mg-dummy", "fe-manual"]
        self._short_name = "solvtype"

class Akeypre(parameter.OneStringParameter):
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
        self._allowed_values = ["unif", "geom"]
        self._short_name = "akeyPRE"

class Akeysolve(parameter.OneStringParameter):
    """ Usage: akeySOLVE {key}

    Specifies how the the finite element mesh should be adaptively subdivided during the
    solve-estimate-refine iterations of a fe-manual finite element calculation. This allows
    for various a posteriori refinement schemes.  key is a text string that specifies the method
    used to guide adaptive refinement:

    * resi - Residual-based a posteriori refinement
    """
    def __init__(self):
        super(Akeysolve, self).__init__()
        self._allowed_values = ["resi"]
        self._short_name = "akeySOLVE"

class Async(parameter.OneIntegerParameter):
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
    def __init__(self):
        super(Async, self).__init__()
        self._short_name = "async"

class Bcfl(parameter.OneStringParameter):
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
        self._allowed_values = ["zero", "sdh", "mdh", "focus", "map"]
        self._short_name = "bcfl"

class Gcent(parameter.Parameter):
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
        super(Gcent, self).__init__()
        for attr in {"mol", "xcent", "ycent", "zcent"}:
            setattr(self, attr, None)
        self._short_name = "gcent"
    def parse(self, tokens):
        token = tokens.pop(0).lower()
        if token == "mol":
            id_token = int(tokens.pop(0))
            self.mol = id_token
        else:
            self.xcent = float(token)
            self.ycent = float(tokens.pop(0))
            self.zcent = float(tokens.pop(0))
    def validate(self):
        # No validation here since it happens in the parsing above
        pass
    def __str__(self):
        outstr = "%s " % self.short_name()
        mol_id = self.mol
        if mol_id:
            outstr = outstr + "mol %d" % mol_id
        else:
            outstr = outstr + "%g %g %g" % (self.xcent, self.ycent, self.zcent)
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
    def __init__(self):
        super(Cgcent, self).__init__()
        self._short_name = "cgcent"

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
    def __init__(self):
        super(Cglen, self).__init__()
        self._short_name = "cglen"

class Chgm(parameter.OneStringParameter):
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
        self._allowed_values = ["spl0", "spl2", "spl4"]
        self._short_name = "chgm"

class Dime(parameter.ThreeIntegerParameter):
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
    def __init__(self):
        super(Dime, self).__init__()
        self._short_name = "dime"
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

class Domainlength(parameter.ThreeFloatParameter):
    """ Specify the rectangular finite element mesh domain lengths for fe-manual finite element
    calculations.  This length may be different in each direction. If the usemesh keyword is
    included, then this command is ignored. The syntax is:

    domainLength {xlen ylen zlen}

    where the parameters xlen, ylen, zlen are floating point numbers that specify the mesh lengths
    in the x-, y-, and z-directions (respectively) in units of &Aring;. """
    def __init__(self):
        super(Domainlength, self).__init__()
        self._short_name = "domainlength"
    def validate(self):
        if (self.xfloat < parameter.FLOAT_EPSILON) \
        or (self.yfloat < parameter.FLOAT_EPSILON) or (self.zfloat < parameter.FLOAT_EPSILON):
            errstr = "One of the domain length parameters is zero or negative"
            raise ValueError(errstr)

class Ekey(parameter.OneStringParameter):
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
        self._allowed_values = ["simp", "global", "frac"]
        self._short_name = "ekey"

class Etol(parameter.OneFloatParameter):
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
    def __init__(self):
        super(Etol, self).__init__()
        self._short_name = "etol"

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
    def __init__(self):
        super(Fgcent, self).__init__()
        self._short_name = "fgcent"

class Fglen(Glen):
    """ Specifies the fine mesh domain lengths in a multigrid focusing calculation (mg-para or
    mg-auto); this may be different in each direction. The syntax is

    fglen {xlen ylen zlen}

    This should enclose the region of interest in the molecule. The arguments to this command are:

    xlen ylen zlen - Grid lengths (floating point numbers) in the x-, y-, and z-directions in
    &Aring;.

    See also: cglen """
    def __init__(self):
        super(Fglen, self).__init__()
        self._short_name = "fglen"

class Ion(parameter.Parameter):
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
        super(Ion, self).__init__()
        self.charge = None
        self.conc = None
        self.radius = None
        self._short_name = "ion"
    def validate(self):
        if self.conc < 0:
            errstr = "Concentration is negative (%g)" % self.conc
            raise ValueError(errstr)
        if self.radius < parameter.FLOAT_EPSILON:
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
    def __str__(self):
        outstr = "ion charge %d conc %g radius %g" % (self.charge, self.conc, self.radius)
        return outstr

class Eqntype(parameter.OneStringParameter):
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
        self._allowed_values = ["lpbe", "lrpbe", "npbe", "nrpbe"]
        self._short_name = "eqntype"

class Maxsolve(parameter.OneIntegerParameter):
    """ Specify the number of times to perform the solve-estimate-refine iteration of the finite
    element solver (fe-manual). The syntax is

    maxsolve { num }

    where num is an integer indicating the desired maximum number of iterations.

    See also: maxvert, targetRes  """
    def __init__(self):
        super(Maxsolve, self).__init__()
        self._short_name = "maxsolve"
    def validate(self):
        if self._parm < 1:
            errstr = "maxsolve is less than 1"
            raise ValueError(errstr)

class Maxvert(parameter.OneIntegerParameter):
    """ Specify the maximum number of vertices to allow during solve-estimate-refine cycle of
    finite element solver (fe-manual). This places a limit on the memory that can be used by the
    solver. The syntax is

    maxvert { num }

    where num is an integer indicating the maximum number of vertices.

    See also: targetRes, maxsolve. """
    def __init__(self):
        super(Maxvert, self).__init__()
        self._short_name = "maxvert"
    def validate(self):
        if self._parm < 1:
            raise ValueError("maxvert is less than 1")

class Nlev(parameter.OneIntegerParameter):
    """ Specify the depth of the multilevel hierarchy used in the mg-manual multigrid solver. See
    dime for a discussion of how nlev relates to grid dimensions. The syntax is

    nlev {lev}

    where lev is an integer indicating the desired depth of the multigrid hierarchy. """
    def __init__(self):
        super(Nlev, self).__init__()
        self._short_name = "nlev"
    def validate(self):
        if self._parm < 1:
            raise ValueError("nlev is less than 1")

class Ofrac(parameter.OneFloatParameter):
    """ Specify the amount of overlap to include between the individual processors meshes in a
    parallel focusing calculation (mg-para). This should be a value between 0 and 1. The syntax is

    ofrac {frac}

    where frac is a floating point value between 0.0 and 1.0 denoting the amount of overlap between
    processors.  Empirical evidence suggests that an ofrac value of 0.1 is sufficient to generate
    stable energies. However, this value may not be sufficient to generate stable forces and/or good
    quality isocontours. """
    def __init__(self):
        super(Ofrac, self).__init__()
        self._short_name = "ofrac"
    def validate(self):
        if self._parm < parameter.FLOAT_EPSILON:
            raise ValueError("ofrac is less than or equal to 0")

class Pdie(parameter.OneFloatParameter):
    """ Specify the dielectric constant of the biomolecule. This is usually a value between 2 to 20,
    where lower values consider only electronic polarization and higher values consider additional
    polarization due to intramolecular motion. The syntax is:

    pdie {diel}

    where diel is the floating point value of the unitless biomolecular dielectric constant.

    See also: sdie """
    def __init__(self):
        super(Pdie, self).__init__()
        self._short_name = "pdie"
    def validate(self):
        if self._parm < 1.0:
            raise ValueError("pdie is less than 1")

class Pdime(parameter.ThreeIntegerParameter):
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
    def __init__(self):
        super(Pdime, self).__init__()
        self._short_name = "pdime"
    def validate(self):
        if (self.xint < 1) or (self.yint < 1) or (self.zint < 1):
            errstr = "One of the dimensions is less than 1"
            raise ValueError(errstr)

class Sdie(parameter.OneFloatParameter):
    """ Specify the dielectric constant of the solvent. Bulk water at biologically-relevant
    temperatures is usually modeled with a dielectric constant of 78-80. The syntax is

    sdie {diel}

    where diel is a floating point number representing the solvent dielectric constant (unitless).

    See also: pdie """
    def __init__(self):
        super(Sdie, self).__init__()
        self._short_name = "sdie"
    def validate(self):
        if self._parm < 1.0:
            raise ValueError("sdie is less than 1.0")

class Targetnum(parameter.OneIntegerParameter):
    """ Specify the target number of vertices in the initial finite element mesh for fe-manual
    calculations.  Initial refinement will continue until this number is reached or the the longest
    edge of every simplex is below targetNum. The syntax is

    targetNum { num }

    where num is an integer denoting the target number of vertices in initial mesh. See also:
    targetRes """
    def __init__(self):
        super(Targetnum, self).__init__(self)
        self._short_name = "targetnum"
    def validate(self):
        if self._parm < 1:
            raise ValueError("targetnum is less than 1")

class Targetres(parameter.OneFloatParameter):
    """ Specify the target resolution of the simplices in a finite element mesh (fe-manual);
    refinement will continue until the longest edge of every simplex is below this value. The syntax
    is

    targetRes { res }

    where res is a floating point number denoting the target resolution for longest edges of
    simplices in mesh (in &Aring;).

    See also: maxvert, maxsolve, targetNum. """
    def __init__(self):
        super(Targetres, self).__init__(self)
        self._short_name = "targetres"
    def validate(self):
        if self._parm < parameter.FLOAT_EPSILON:
            raise ValueError("targetres is zero or negative")

class Usemap(parameter.Parameter):
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
        self._allowed_values = ["diel", "kappa", "charge", "pot"]
        self.map_id = None
        self._short_name = "usemap"
    def parse(self, tokens):
        token = tokens.pop(0).lower()
        if token in self._allowed_values:
            self.type = token
        else:
            errstr = "Unknown map type (%s) in usemap" % token
            raise ValueError(errstr)
        self.map_id = int(tokens.pop(0))
    def validate(self):
        # Validation happens in parsing
        pass
    def __str__(self):
        outstr = "usemap %s %d" % (self.type, self.map_id)
        return outstr

class Usemesh(parameter.OneIntegerParameter):
    """ Specify the external finite element mesh to be used in the finite element Poisson-Boltzmann
    calculation (fe-manual). These must have been input via an earlier READ mesh statement. The
    syntax is

    usemesh {id}

    where id is an integer ID specifying the particular map read in with READ mesh. These IDs are
    assigned sequentially, starting from 1, and incremented independently for each mesh read by
    APBS. """
    def __init__(self):
        super(Usemesh, self).__init__()
        self._short_name = "usemesh"

class Write(parameter.Parameter):
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
        self._short_name = "write"
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
    def __str__(self):
        return "write %s %s \"%s\"" % (self.type, self.format, self.stem)

class Writemat(parameter.FormatPathParameter):
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
        self._allowed_values = ["poisson"]
        self._short_name = "writemat"
    def parse(self, tokens):
        self.type = tokens.pop(0).lower()
        self.stem = tokens.pop(0)
    def validate(self):
        if not self.type in self._allowed_values:
            errstr = "Unknown token (%s) for %s" % self.type
            raise ValueError(errstr)
    def __str__(self):
        return "writemat %s %s\n" % (self.type, self.stem)

class Elec(parameter.ParameterSection):
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
        self._short_name = "elec"
        self._allowed_keywords = {"name" : Name, "solvtype" : Solvtype, "akeypre" : Akeypre,
                                 "akeysolve" : Akeysolve, "async" : Async, "bcfl" : Bcfl,
                                 "calcenergy" : Calcenergy, "calcforce" : Calcforce,
                                 "cgcent" : Cgcent, "cglen" : Cglen, "chgm" : Chgm, "dime" : Dime,
                                 "domainlength" : Domainlength, "ekey" : Ekey, "etol" : Etol,
                                 "fgcent" : Fgcent, "fglen" : Fglen, "gcent" : Gcent, "glen" : Glen,
                                 "grid" : Grid, "eqntype" : Eqntype, "maxsolve" : Maxsolve,
                                 "maxvert" : Maxvert, "mol" : Mol, "nlev" : Nlev, "ofrac" : Ofrac,
                                 "pdie" : Pdie, "pdime" : Pdime, "sdens" : Sdens, "sdie" : Sdie,
                                 "srad" : Srad, "srfm" : Srfm, "swin" : Swin,
                                 "targetnum" : Targetnum, "targetres" : Targetres, "temp" : Temp,
                                 "usemap" : Usemap, "usemesh" : Usemesh, "write" : Write,
                                 "writemat" : Writemat}
        # Ions are different because you can have many of them
        self._allowed_keywords["ion"] = None
    def parse(self, tokens):
        token = tokens.pop(0)
        while True:
            token_name = token.lower()
            _LOGGER.debug(token_name)
            if token_name in ["ion"]:
                ion = Ion()
                ion.parse(tokens)
                if hasattr(self, "ion"):
                    self.ion.append(ion)
                else:
                    self.ion = [ion]
            elif token_name in Solvtype()._allowed_values:
                # This is a special section to handle the old-format ELEC solver type declaration
                solvtype = Solvtype()
                solvtype.parse([token_name])
                solvtype.validate()
                self.solvtype = solvtype
                _LOGGER.debug(solvtype)
            elif token_name in Eqntype()._allowed_values:
                # This is a special section to handle the old-format ELEC equation type declaration
                eqntype = Eqntype()
                eqntype.parse([token_name])
                eqntype.validate()
                self.eqntype = eqntype
                _LOGGER.debug(eqntype)
            elif token_name == "end":
                return
            elif token_name in self._allowed_keywords:
                    token_object = self._allowed_keywords[token_name]()
                    token_object.parse(tokens)
                    token_object.validate()
                    setattr(self, token_name, token_object)
            else:
                errstr = "Unknown token (%s) for %s" % (token, self.short_name())
                raise ValueError(errstr)
            token = tokens.pop(0)
    def validate_mgauto(self):
        """ Validate current input file contents """
        required_parms = ["bcfl", "chgm", "cgcent", "cglen", "dime", "fgcent", "fglen",
                          "eqntype", "mol", "pdie", "sdens", "sdie", "srad", "srfm", "swin",
                          "temp"]
        for parm in required_parms:
            if not parm in self.contents():
                errstr = "Missing required parameter %s for %s" % (parm, self.solvtype)
                raise ValueError(errstr)
    def validate_mgpara(self):
        """ Validate current input file contents """
        required_parms = ["bcfl", "chgm", "cgcent", "cglen", "dime", "fgcent", "fglen",
                          "eqntype", "mol", "ofrac", "pdie", "pdime", "sdens", "sdie", "srad",
                          "srfm", "swin", "temp"]
        for parm in required_parms:
            if not parm in self.contents():
                errstr = "Missing required parameter %s for %s" % (parm, self.solvtype)
                raise ValueError(errstr)
    def validate_mgmanual(self):
        """ Validate current input file contents """
        required_parms = ["bcfl", "chgm", "dime", "gcent", "eqntype", "mol", "nlev", "pdie",
                          "sdens", "sdie", "srad", "srfm", "swin", "temp"]
        for parm in required_parms:
            if not parm in self.contents():
                errstr = "Missing required parameter %s for %s" % (parm, self.solvtype)
                raise ValueError(errstr)
        if (not "glen" in self.contents()) and (not "grid" in self.contents()):
            errstr = "Either grid or glen needs to be specified for mg-manual"
            raise ValueError(errstr)
        if self.solvtype == "mg-manual":
            dime = self.dime
            nlev = self.nlev
            ncurr = nlev.parm()
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
            if not parm in self.contents():
                errstr = "Missing required parameter %s for %s" % (parm, self.solvtype)
                raise ValueError(errstr)
    def validate(self):
        """ Validate current input file contents """
        solvtype = self.solvtype.parm()
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
        keys = self.contents()
        keys.remove("name")
        outstr = outstr + "\t%s\n" % self.solvtype
        keys.remove("solvtype")
        outstr = outstr + "\t%s\n" % self.eqntype
        keys.remove("eqntype")
        for key in keys:
            values = getattr(self, key)
            try:
                for value in values:
                    outstr = outstr + "\t%s\n" % value
            except TypeError:
                outstr = outstr + "\t%s\n" % values
        outstr = outstr + "end\n"
        return outstr
