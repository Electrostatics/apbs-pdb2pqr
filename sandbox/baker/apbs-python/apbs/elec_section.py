""" Parse APBS input file ELEC sections """
from parameters import *

class Name(OneStringParameter):
    """ Usage: name {id}
        
    Since numerous ELEC blocks may appear in an APBS input file, it can be difficult to keep
    track of them all. It is possible to assign an optional name to each ELEC block to simplify
    the organizational process.
        
    {id} is an alphanumeric string denoting the "name" of the calculation block."""
    @property
    def name(self):
        return "name"

class Type(OneStringParameter):
    """ Usage:  type {type}
        
    Specify the type of electrostatics calculation to be performed where {type} is one of:
        
    * fe-manual - manual finite element method
    * mg-auto - automatic multigrid setup
    * mg-dummy - dummy multigrid setup for writing out maps
    * mg-manual - manual multigrid setup
    * mg-para - parallel multigrid setup
        
    This isn't actually official yet.  I would like to change the default behavior of APBS
    from stating the ELEC calculation type without any parameters to adding the "type" keyword."""
    allowed_values = ["mg-auto", "mg-para", "mg-manual", "fe-manual", "mg-dummy"]
    @property
    def name(self):
        return "type"

class Akeypre(OneStringParameter):
    """ Usage: akeyPRE {key}
        
    Specifies how the initial finite element mesh should be constructed (from refinement of a
    very coarse 8-tetrahedron mesh prior to the solve-estimate-refine iteration in fe-manual
    finite element calculations.
        
    key is a text string that specifies the method used to guide initial refinement and takes
    one of the values:
        
    * unif - Uniform refinement 
    * geom - Geometry-based refinement at molecular surfaces and charges """
    allowed_values = ["unif", "geom"]
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
    allowed_values = ["resi"]
    @property
    def name(self):
        return "akeysolve"

class Elec(Parameter):
    """ ELEC section for APBS input file """
    def __init__(self):
        self.content_dict = {}
        
    @property
    def name(self):
        return "elec"
    
    def parse(self, tokens):
        token = tokens.pop(0)
        while True:
            tokenName = token.lower()
            if tokenName == "name":
                name = Name()
                name.parse(tokens)
                name.validate()
                self.content_dict["name"] = name
            elif tokenName == "type":
                elecType = Type()
                elecType.parse(tokens)
                elecType.validate()
                self.content_dict["type"] = elecType
            elif tokenName in Type.allowed_values:
                elecType = Type()
                elecType.parse([tokenName])
                elecType.validate()
                self.content_dict["type"] = elecType
            elif tokenName == "akeypre":
                akeypre = Akeypre()
                akeypre.parse(tokens)
                akeypre.validate()
                self.content_dict["akeypre"] = akeypre
            elif tokenName == "akeysolve":
                akeysolve = Akeysolve()
                akeysolve.parse(tokens)
                akeysolve.validate()
                self.content_dict["akeysolve"] = akeysolve
            else:
                errstr = "Unknown token (%s) for %s" % (token, self.name)
                raise ValueError, errstr
            token = tokens.pop(0)
    def validate(self):
        import sys
        sys.stderr("ELEC validation incomplete!\n")
        
    def contents(self):
        return self.content_dict
        
    def __str__(self):
        outstr = "elec\n"
        keys = self.content_dict.keys()
        try:
            outstr = outstr + "\t%s\n" % self.content_dict["name"]
            keys.remove("name")
        except KeyError:
            pass
        outstr = outstr + "\t%s\n" % self.content_dict["type"]
        keys.remove("type")
        for key in keys:
            outstr = outstr + "\t%s\n" % self.content_dict[key]
        outstr = outstr + "end\n"
        return outstr
    

    
    def parse_async(self):
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
        
        rankToken = self.tokens.pop(0)
        try:
            rank = int(rankToken)
            self.parameters["async"] = rank
        except ValueError:
            parameterError(tokenName, rankToken)
    
    def parse_bcfl(self):
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
        typeName = self.tokens.pop(0)
        typeName = typeName.lower()
        if typeName in ["zero", "sdh", "mdh", "focus", "map"]:
            self.parameters["bcfl"] = typeName
        else:
            parameterError(tokenName, typeName)
    
    def parse_calcenergy(self):
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
        typeName = self.tokens.pop(0)
        typeName = typeName.lower()
        # TODO - This keyword is stupid; why don't we just calculate this by default?
        if typeName in ["no", "total", "comps"]:
            self.parameters["calcenergy"] = typeName
        else:
            parameterError(tokenName, typeName)
    
    def parse_calcforce(self):
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
        typeName = self.tokens.pop(0)
        typeName = typeName.lower()
        # TODO - This keyword is stupid; why don't we just calculate this by default?
        if typeName in ["no", "total", "comps"]:
            self.parameters["calcforce"] = typeName
        else:
            parameterError(tokenName, typeName)
    
    def parseCRAP(self):
        """ Parse input file """
        token = self.tokens.pop(0)
        while True:
            tokenName = token.lower()
            if tokenName == "elec":
                pass
            elif tokenName == "name":
                self.parameters["name"] = Name(self.tokens)
            elif tokenName == "type":
                self.parse_type()
            elif tokenName in ["mg-auto", "mg-para", "mg-manual", "fe-manual", "mg-dummy"]:
                self.parse_type(tokenName)
            elif tokenName == "akeypre":
                self.parse_akeypre()
            elif tokenName == "akeysolve":
                self.parse_akeysolve()
            elif tokenName == "async":
                self.parse_async()
            elif tokenName == "bcfl":
                self.parse_bcfl()
            elif tokenName == "calcenergy":
                self.parse_calcenergy()
            elif tokenName == "calcforce":
                self.parse_calcforce()
            elif tokenName == "cgcent":
                self.parse_cgcent()
            else:
                errstr = "Unknown ELEC token (%s)" % tokenName
                raise ValueError, errstr
            token = self.tokens.pop(0)
        