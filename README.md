APBS and PDB2PQR: Electrostatic and solvation properties from complex molecules.
============

Please visit our home page for more information: [www.poissonboltzmann.org](http://www.poissonboltzmann.org).

## Quick links
* Binary releases can be found at [GitHub](https://github.com/Electrostatics/apbs-pdb2pqr/releases) and [SourceForge](https://sourceforge.net/projects/apbs/). 
* Source code build instructions are available for [APBS](https://github.com/Electrostatics/apbs-pdb2pqr/blob/master/apbs/README.md) and [PDB2PQR](https://github.com/Electrostatics/apbs-pdb2pqr/blob/master/pdb2pqr/README.md)
* Support can be obtained through [mailing lists](http://www.poissonboltzmann.org/support/home/) or [chat](https://gitter.im/Electrostatics/help)
* Contributing software we include as git submodules:
  * [PB-AM](https://github.com/davas301/pb_solvers)
  * [PB-SAM](https://github.com/davas301/pb_solvers)
  * [TABI](https://github.com/lwwilson1/TABIPB)
  * [GeoFlow](https://github.com/Electrostatics/geoflow_c)

## More information about the software

An understanding of electrostatic interactions is essential for the study of biomolecular processes. The structures of proteins and other biopolymers are being determined at an increasing rate through structural genomics and other efforts while specific linkages of these biopolymers in cellular pathways or supramolecular assemblages are being detected by genetic and proteomic studies. To integrate this information in physical models for drug discovery or other applications requires the ability to evaluate the energetic interactions within and between biopolymers. Among the various components of molecular energetics, solvation properties and electrostatic interactions are of special importance due to the long range of these interactions and the substantial charges of typical biopolymer components. 

APBS is a unique software package which solves the equations of continuum electrostatics for large biomolecular assemblages. This software was designed “from the ground up” using modern design principles to ensure its ability to interface with other computational packages and evolve as methods and applications change over time. The APBS code is accompanied by extensive documentation for both users and programmers and is supported by a variety of utilities for preparing calculations and analyzing results. Finally, the free, open-source APBS license ensures its accessibility to the entire biomedical community. 

The use of continuum solvation methods such as APBS requires accurate and complete structural data as well as force field parameters such as atomic charges and radii. Unfortunately, the limiting step in continuum electrostatics calculations is often the addition of missing atomic coordinates to molecular structures from the Protein Data Bank and the assignment of parameters to these structures. To address this problem, we have developed PDB2PQR. This software automates many of the common tasks of preparing structures for continuum solvation calculations as well as many other types of biomolecular structure modeling, analysis, and simulation. These tasks include: adding a limited number of missing heavy (non-hydrogen) atoms to biomolecular structures, estimating titration states and protonating biomolecules in a manner consistent with favorable hydrogen bonding, assigning charge and radius parameters from a variety of force fields, and finally generating “PQR” output compatible with several popular computational biology modeling and analysis packages. This service is intended to facilitate the setup and execution of electrostatics calculations for both experts and non-experts and broaden the accessibility of biomolecular solvation/electrostatics analyses to the research community. 
