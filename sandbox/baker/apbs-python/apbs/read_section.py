""" Parse the READ input file section """
from parameters import *

# Set the following flag to True to use deprecated input format
useDeprecated = True

class Read:
    """ READ input file section """
    def __init__(self):
        self.content_dict = { "diel" : [], "kappa" : [], "mesh" : [], "mol" : [], "parm" : [], "pot" : []}
    
    @property
    def name(self):
        return "read"

    def validate(self):
        pass
    
    def contents(self):
        return self.content_dict
    
    def __str__(self):
        outstr = "read\n"
        for key, values in self.content_dict.items():
            for value in values:
                outstr = outstr + "\t%s\n" % value
        outstr = outstr + "end\n"
        return outstr

    def parse(self, tokens):
        """ Parse input file """
        token = tokens.pop(0)
        while True:
            tokenName = token.lower()
            if tokenName == "charge":
                charge = Charge()
                charge.parse(tokens)
                charge.validate()
                self.content_dict["charge"].append(charge)
            elif tokenName == "diel":
                diel = Diel()
                diel.parse(tokens)
                diel.validate()
                self.content_dict["diel"].append(diel)
            elif tokenName == "kappa":
                kappa = Kappa()
                kappa.parse(tokens)
                kappa.validate()
                self.content_dict["kappa"].append(kappa)
            elif tokenName == "mesh":
                mesh = Mesh()
                mesh.parse(tokens)
                mesh.validate()
                self.content_dict["mesh"].append(mesh)
            elif tokenName == "mol":
                mol = Mol()
                mol.parse(tokens)
                mol.validate()
                self.content_dict["mol"].append(mol)
            elif tokenName == "parm":
                parm = Parm()
                parm.parse(tokens)
                parm.validate()
                self.content_dict["parm"].append(parm)
            elif tokenName == "pot":
                pot = Pot()
                pot.parse(tokens)
                pot.validate()
                self.content_dict["pot"].append(pot)
            elif tokenName == "end":
                return
            else:
                errstr = "Unknown READ token (%s)" % tokenName
                raise ValueError, errstr
            token = tokens.pop(0)

class Charge(Parameter):
    """ Charge input file
    
    Usage:  charge {format} {path}
        
        This command allows APBS to read the fixed (molecular) charge density function mapped to a mesh.
        The inputs are maps of charge densities; these values have units of ec &Aring;^-3, where ec is the
        electron charge. In general, this command will read charge-maps written by write commands in
        earlier APBS calculations.
        
        Arguments for this command are:
        
        * format - The format of the charge map. Acceptable values include:
          - dx - OpenDX format 
          - gz - gzipped (zlib) compressed OpenDX format. Files can be read directly in compressed form.
        * path - The location of the charge map file.
    """
    def __init__(self):
        self.content_dict = { "format" : None, "path" : None}
    
    @property
    def name(self):
        return "charge"
    
    def contents(self):
        return self.content_dict

    def parse(self, tokens):
        self.content_dict["format"] = tokens.pop(0)
        self.content_dict["path"] = tokens.pop(0)
    
    def validate(self):
        format = self.content_dict["format"].lower()
        if not format in ["dx", "dz"]:
            errstr = "Unknown format %s for %s" % (format, self.name)
            raise ValueError, errstr
        if not self.content_dict["path"]:
            errstr = "Missing path for %s" % self.name
            raise ValueError, errstr
    
    def __str__(self):
        return " ".join([self.name, self.content_dict["format"], self.content_dict["path"]])

class Diel(Parameter):
    """ Dielectric input file
    
    Usage:  diel {format} {path-x} {path-y} {path-z}
        
        This command allows APBS to read the dielectric function mapped to 3 meshes shifted by
        one-half grid spacing in the x, y, and z directions. The inputs are maps of dielectric
        variables between the solvent and biomolecular dielectric constants; these values are
        unitless. In general, this command will read dielectric maps written by write commands in
        earlier APBS calculations.
        
        NOTE: if you choose this option and have a non-zero ionic strength, you must also include a read kappa statement
        
        Required arguments for this command are:
        
        * format - The format of the dielectric map. Acceptable values include:
          - dx - OpenDX format
          - gz - gzipped (zlib) compressed OpenDX format. Files can be read directly in compressed form.
        * path-x - The location of the x-shifted dielectric map file. 
        * path-y - The location of the y-shifted dielectric map file. 
        * path-z - The location of the z-shifted dielectric map file. """
    def __init__(self):
        self.content_dict = {"format" : None, "path-x" : None, "path-y" : None, "path-z" : None }
    
    @property
    def name(self):
        return "diel"
    
    def contents(self):
        return self.content_dict

    def parse(self, tokens):
        self.content_dict["format"] = tokens.pop(0)
        self.content_dict["path-x"] = tokens.pop(0)
        self.content_dict["path-y"] = tokens.pop(0)
        self.content_dict["path-z"] = tokens.pop(0)
    
    def validate(self):
        format = self.content_dict["format"].lower()
        if not format in ["dx", "gz"]:
            errstr = "Unknown format %s for %s" % (format, self.name)
            raise ValueError, errstr
        if not (self.content_dict["path-x"] or self.content_dict["path-y"] or self.content_dict["path-y"]):
            errstr = "Missing input path for %s" % self.name
            raise ValueError, errstr
    
    def __str__(self):
        return " ".join([self.name, self.content_dict["format"], self.content_dict["path-x"],
                         self.content_dict["path-y"], self.content_dict["path-z"]])

class Kappa(Parameter):
    """ Kappa input file
    
    Usage:  kappa {format} {path}
    
    This command allows APBS to read the ion-accessibility function mapped to a mesh. The
    inputs are maps of ion accessibility values which range between 0 and the build
    Debye-Huckel screening parameter; these values have units of &Aring;^-2. In general, this
    command will read kappa-maps written by write commands in earlier APBS calculations.
    
    NOTE: if you choose this option, you must also include a read diel statement. 
    
    Arguments for this command are:

    * format - The format of the kappa map. Acceptable values include:
      - dx - OpenDX format
      - gz - gzipped (zlib) compressed OpenDX format. Files can be read directly in compressed form.
    * path - The location of the kappa map file."""
    def __init__(self):
        self.content_dict = { "format" : None, "path" : None}
    
    @property
    def name(self):
        return "kappa"
    
    def parse(self, tokens):
        self.content_dict["format"] = tokens.pop(0)
        self.content_dict["path"] = tokens.pop(0)
    
    def contents(self):
        return self.content_dict

    def validate(self):
        format = self.content_dict["format"].lower()
        if not format in ["dx", "gz"]:
            errstr = "Unknown format %s for %s" % (format, self.name)
            raise ValueError, errstr
        if not (self.content_dict["path"]):
            errstr = "Missing input path for %s" % self.name
            raise ValueError, errstr
    
    def __str__(self):
        return " ".join([self.name, self.content_dict["format"], self.content_dict["path"]])

class Mesh(Parameter):
    """ Mesh input file
    
    Usage:  mesh {format} {path}
        
    This command allows APBS to read a finite element mesh to use as a starting point for finite
    element calculations. The input is simply the mesh geometry; e.g., as produced by a finite
    element mesh generation program such as LBIE-Mesher or GAMer. Arguments for this command are:
    
    * format - The format of the input mesh. Acceptable values include:
      - mcsf - MCSF format 
    * path - The location of the mesh file. """
    def __init__(self):
        self.content_dict = { "format" : None, "path" : None}
    
    @property
    def name(self):
        return "mesh"
    
    def parse(self, tokens):
        self.content_dict["format"] = tokens.pop(0)
        self.content_dict["path"] = tokens.pop(0)

    def contents(self):
        return self.content_dict
    
    def validate(self):
        format = self.content_dict["format"].lower()
        if not format in ["mcsf"]:
            errstr = "Unknown format %s for %s" % (format, self.name)
            raise ValueError, errstr
        if not (self.content_dict["path"]):
            errstr = "Missing input path for %s" % self.name
            raise ValueError, errstr
    
    def __str__(self):
        return " ".join([self.name, self.content_dict["format"], self.content_dict["path"]])

class Mol(Parameter):
    """ Information about APBS molecule files
    
    Usage:  mol {format} {path} 
    
    This command specifies the molecular data to be read into APBS. The required arguments are:
    
    * format - The format of the input data. Acceptable values include:
      - pqr - Specify that molecular data is in PQR format
      - pdb - Specify that molecular data is in pseudo-PDB format.  If this type of structure
      file is used, then a parameter file must also be specified to provide charge and radius
      parameters for the biomolecule's atoms. 
    * path - The location of the molecular data file. """
    def __init__(self):
        self.content_dict = { "format" : None, "path" : None}
    
    @property
    def name(self):
        return "mol"

    def contents(self):
        return self.content_dict
    
    def parse(self, tokens):
        self.content_dict["format"] = tokens.pop(0)
        self.content_dict["path"] = tokens.pop(0)
    
    def validate(self):
        format = self.content_dict["format"].lower()
        if not format in ["pdb", "pqr"]:
            errstr = "Unknown format %s for %s" % (format, self.name)
            raise ValueError, errstr
        if not (self.content_dict["path"]):
            errstr = "Missing input path for %s" % self.name
            raise ValueError, errstr
    
    def __str__(self):
        return " ".join([self.name, self.content_dict["format"], self.content_dict["path"]])

class Parm(Parameter):
    """ Information about APBS parameter files
    
    Usage: parm {format} {path}
    
    This command specifies the charge and radius data to be used with pseudo-PDB-format molecule
    files. The arguments are:
    
    * format - The format of the parameter file. Acceptable flags include:
      - flat - Specify that the parameter file is in APBS Flat-file parameter format 
      - xml - Specify that the parameter file is in APBS XML parameter format 
    * path - The location of the parameter data file. 

    Note that APBS provides a few example files as part of the source code distribution.
    Currently, example files only contain the polar parameters that can also be assigned more
    easily through the PDB2PQR software.  Parameter files with apolar values are not currently
    available for protein and nucleic acid parameters and are actively under development as a
    research project.  Please contact Nathan Baker for additional information about the state of
    this research, particularly if you are interested in helping. """
    def __init__(self):
        self.content_dict = { "format" : None, "path" : None}
    
    @property
    def name(self):
        return "parm"

    def contents(self):
        return self.content_dict

    def parse(self, tokens):
        self.content_dict["format"] = tokens.pop(0)
        self.content_dict["path"] = tokens.pop(0)
    
    def validate(self):
        format = self.content_dict["format"].lower()
        if not format in ["flat", "xml"]:
            errstr = "Unknown format %s for %s" % (format, self.name)
            raise ValueError, errstr
        if not (self.content_dict["path"]):
            errstr = "Missing input path for %s" % self.name
            raise ValueError, errstr
    
    def __str__(self):
        return " ".join([self.name, self.content_dict["format"], self.content_dict["path"]])

class Pot(Parameter):
    """ Information about APBS potential files
    
    Usage:  pot {format} {path}
        
    This command allows APBS to read the electrostatic potential mapped to a mesh. The inputs
    are maps of the electrostatic potential from a previous calculation. In general, this command
    will read potential-maps written by write commands in earlier APBS calculations.
    
    NOTE: To use this functionality you must set the bcfl keyword to map. See also: usemap
    
    Arguments for this command are:
    * format - The format of the potential map. Acceptable values include:
      - dx - OpenDX format
      - gz - gzipped (zlib) compressed OpenDX format. Files can be read directly in compressed form.
    * path - The location of the potential map file."""
    def __init__(self):
        self.content_dict = { "format" : None, "path" : None}
    
    @property
    def name(self):
        return "pot"

    def contents(self):
        return self.content_dict
    
    def parse(self, tokens):
        self.content_dict["format"] = tokens.pop(0)
        self.content_dict["path"] = tokens.pop(0)
    
    def validate(self):
        format = self.content_dict["format"].lower()
        if not format in ["dx", "gz"]:
            errstr = "Unknown format %s for %s" % (format, self.name)
            raise ValueError, errstr
        if not (self.content_dict["path"]):
            errstr = "Missing input path for %s" % self.name
            raise ValueError, errstr
    
    def __str__(self):
        return " ".join([self.name, self.content_dict["format"], self.content_dict["path"]])

