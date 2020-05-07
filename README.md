<a href='http://www.poissonboltzmann.org/'><img src='https://img.shields.io/badge/homepage-poissonboltzmann-blue.svg'></a>

[![Build Status](https://travis-ci.org/Electrostatics/apbs-pdb2pqr.svg?branch=master)](https://travis-ci.org/Electrostatics/apbs-pdb2pqr)

# APBS and PDB2PQR: electrostatic and solvation properties for complex molecules

## Getting started

* Please [register your use](http://eepurl.com/by4eQr) to ensure continued support for the APBS-PDB2PQR software.
* Read our [online documentation](http://apbs-pdb2pqr.readthedocs.io/).
* Get started with the [web server](http://server.poissonboltzmann.org/).
* Download the software [following these instructions](http://apbs-pdb2pqr.readthedocs.io/en/latest/downloads.html).

-----


## More information about the software

An understanding of electrostatic interactions is essential for the study of biomolecular processes.
The structures of proteins and other biopolymers are being determined at an increasing rate through structural genomics and other efforts while specific linkages of these biopolymers in cellular pathways or supramolecular assemblages are being detected by genetic and proteomic studies.
To integrate this information in physical models for drug discovery or other applications requires the ability to evaluate the energetic interactions within and between biopolymers.
Among the various components of molecular energetics, solvation properties and electrostatic interactions are of special importance due to the long range of these interactions and the substantial charges of typical biopolymer components.

### APBS

APBS solves the equations of continuum electrostatics for large biomolecular assemblages.
This software was designed “from the ground up” using modern design principles to ensure its ability to interface with other computational packages and evolve as methods and applications change over time.
The APBS code is accompanied by extensive documentation for both users and programmers and is supported by a variety of utilities for preparing calculations and analyzing results.
Finally, the free, open-source APBS license ensures its accessibility to the entire biomedical community.

### PDB2PQR
The use of continuum solvation methods such as APBS requires accurate and complete structural data as well as force field parameters such as atomic charges and radii.
Unfortunately, the limiting step in continuum electrostatics calculations is often the addition of missing atomic coordinates to molecular structures from the [Protein Data Bank](http://www.wwpdb.org/) and the assignment of parameters to these structures.
To address this problem, we have developed PDB2PQR.
This software automates many of the common tasks of preparing structures for continuum solvation calculations as well as many other types of biomolecular structure modeling, analysis, and simulation.
These tasks include:

* adding a limited number of missing heavy (non-hydrogen) atoms to biomolecular structures,
* estimating titration states and protonating biomolecules in a manner consistent with favorable hydrogen bonding,
* assigning charge and radius parameters from a variety of force fields, and
* generating “PQR” output compatible with several popular computational biology modeling and analysis packages.

This service is intended to facilitate the setup and execution of electrostatics calculations for both experts and non-experts and broaden the accessibility of biomolecular solvation/electrostatics analyses to the research community.

### Support for APBS-PDB2PQR

APBS and PDB2PQR are supported by NIH grant GM69702.
Additional support and contributors are listed in the [online documentation](http://apbs-pdb2pqr.readthedocs.io/).

### APBS Datasheet

OS | PYTHON VERSION | GEOFLOW | BEM,MSMS | FETK | PBSAM | PBAM | PYTHON | SHARED_LIBS | TESTS PASS
------------- | ------------ | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- 
Ubuntu latest | 3.6, 3.7 | :grinning: | :grinning: | :grinning: | :nauseated_face: | :nauseated_face: | :grinning: | :grinning: | :partying_face:
MacOSX latest | 3.6, 3.7 | :grinning: | :grinning: | :grinning: | :nauseated_face: | :nauseated_face: | :grinning: | :grinning: | :partying_face: