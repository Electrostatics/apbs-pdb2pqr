[![Homepage](https://img.shields.io/badge/homepage-poissonboltzmann-blue.svg)](http://www.poissonboltzmann.org)
[![Docs Build](https://readthedocs.org/projects/apbs-pdb2pqr/badge/)](https://apbs-pdb2pqr.readthedocs.io/)
[![Travis Build](https://travis-ci.org/Electrostatics/apbs-pdb2pqr.svg?branch=master)](https://travis-ci.org/Electrostatics/apbs-pdb2pqr)
![GitHub Build](https://github.com/Electrostatics/apbs-pdb2pqr/workflows/Build/badge.svg)
![Appveyor Build](https://ci.appveyor.com/api/projects/status/github/Electrostatics/apbs-pdb2pqr?branch=master&svg=true)

# APBS: electrostatic and solvation properties for complex molecules

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

APBS solves the equations of continuum electrostatics for large biomolecular assemblages.
This software was designed “from the ground up” using modern design principles to ensure its ability to interface with other computational packages and evolve as methods and applications change over time.
The APBS code is accompanied by extensive documentation for both users and programmers and is supported by a variety of utilities for preparing calculations and analyzing results.
Finally, the free, open-source APBS license ensures its accessibility to the entire biomedical community.

### Support for APBS

APBS is supported by NIH grant GM69702.
Additional support and contributors are listed in the [online documentation](http://apbs-pdb2pqr.readthedocs.io/).

### Platform support

This shows the status of APBS functionality on different platforms.

OS            | PYTHON VERSION | GEOFLOW  | BEM, MSMS | FETK  | PBSAM | PBAM | PYTHON | SHARED_LIBS |
------------- | -------------- | -------- | --------- | ----- | ----- | ---- | ------ | ----------- |
Ubuntu latest | 3.6, 3.7       | ✔️       | ✔️        | ✔️   | ✔️    | ✔️  | ✔️    | ✔️          |
MacOSX latest | 3.6, 3.7       | ✔️       | ✔️        | ✔️   | ✔️    | ✔️  | ✔️    | ✔️          |
Windows 10    | 3.7            | ✔️       | ✔️        | ❌   | ❌    | ✔️  | ✔️    | ❌          |
