#Build Configuration file for pdb2pqr
#While you can specify these on the command line with var=value
# this is the recommended way to setup a build.

#Uncomment the values you would like to change and set new values.


#Installation PREFIX
#Sets the install location of pdb2pqr.
#This defaults to ~/pdb2pqr

#PREFIX="~/pdb2pqr"


#Website URL
#Sets the url of pdb2pqr's web interface.
#This defaults to http://<COMPUTER NAME>/pdb2pqr/
# where <COMPUTER NAME> is the network name of the computer.

#URL="http://<COMPUTER NAME>/pdb2pqr/"


#APBS binary
#Change this to specify the location of the APBS binary if installed.

#APBS=""


#OPAL service URL
#Set this value to use an opal service for processing 
# pdb2pqr jobs for the web front end.

#OPAL="http://nbcr-222.ucsd.edu/opal2/services/pdb2pqr_1.8"


#APBS_OPAL service URL
#Set this value to use an opal service for processing 
# apbs jobs for the web front end.

#APBS_OPAL="http://nbcr-222.ucsd.edu/opal2/services/apbs_1.3"


#MAX_ATOMS
#Sets the maximum number of atoms in a protein for non-Opal job submission.
#Only affects web tools. Default is 10000.

#MAX_ATOMS=10000


#BUILD_PDB2PKA
#Set to False to skip building ligand and pdb2pka support. Requires numpy.
# Defaults to True

#BUILD_PDB2PKA=False


#REBUILD_SWIG
#Set to True to rebuild the swig bindings. 
# Requires swig on the the user path.
# Defaults to False

#REBUILD_SWIG=True