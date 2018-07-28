APBS-PDB2PQR overview
=====================

=====================
Why APBS and PDB2PQR?
=====================

An understating of electrostatic interactions is essential for the study of bio-molecular processes.
The structures of proteins and other bio-polymers are being determined at an increasing rate through structural genomics and other efforts while specific linkages of these biopolymers in cellular pathways or supramolecular assemblages are being detected by genetic and proteomic studies.
To integrate this information in physical models for drug discovery or other applications requires the ability to evaluate the energetic interactions within and between bio-polymers.
Among the various components of molecular energetics, solvation properties and electrostatic interactions are of special importance due to the long range of these interactions and the substantial charges of typical biopolymer components.

APBS is a unique software which solves the equations of continuum electrostatics for large biomolecular assemblages. This software was designed "from the ground up" using modern design principles to ensure its ability to interface with other computational packages and evolve as methods and applications change over time. The APBS code is accompanied by extensive documentation for both users and programmers and is supported by a variety of utilities for preparing calculations and analyzing results. Finally, the free, open-source APBS license ensures its accessibility to the entire biomedical community.

The use of continuum solvation methods such as APBS requires accurate and complete structural data as well as force field parameters such as atomic charges and radii. Unfortunately, the limiting step in continuum electrostatics calculations is often the addition of missing atomic coordinates to molecular structures from the Protein Data Bank and the assignment of parameters to these structures. To adds this problem, we have developed PDB2PQR. This software automates many of the common tasks of preparing structures for continuum solvation calculations as well as many other types of biomolecular structure modeling, analysis, and simulation. These tasks include:

* Adding a limited number of missing heavy (non-hydrogen) atoms to biomolecular structures.
* Estimating titration states and protonating biomolecules in a manner consistent with favorable hydrogen bonding.
* Assigning charge and radius parameters from a variety of force fields.
* Generating "PQR" output compatible with several popular computational modeling and analysis packages.

This service is intended to facilitate the setup and execution of electrostatics calculations for both experts and non-experts and thereby broaden the accessibility of biomolecular solvation and electrostatics analyses to the biomedical community.

Both APBS and PDB2PQR have enjoyed widespread adoption throughout the biomedical community and are used in numerous applications involving biomolecular structures.

==================
What do I do next?
==================

#. **Before you start**, please make sure to support the APBS/PDB2PQR team by `registering your use of the software <http://eepurl.com/by4eQr>`_.
#. **Download the software** following the instructions at :doc:`downloads`.

If you come across anything along the way that we haven’t covered, is incorrect information, or if you know of a tip you think others would find useful, please let us know (see :doc:`getting-help`) and we’ll work on including it in this guide.
