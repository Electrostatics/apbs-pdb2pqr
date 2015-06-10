apbs-pdb2pqr
============

This is the home for the [APBS and PDB2PQR software](http://www.poissonboltzmann.org) source repository.  Binary downloads can be found on [SourceForge](https://sourceforge.net/projects/apbs/). 

# Latest release

If you are looking for the source for the latest release of APBS, be sure to pull the 1.4.1 release branch:

`git checkout 1.4.1-binary-release`

# Geometric Flow

If you want to use the geometric flow implementation, do the following:

1. Get the submodule from github:
  * git submodule init
  * git submodule update
2. In CMake:
  * Set ENABLE_GEOFLOW to ON
3. Build as usual
