# APBS and PDB2PQR: electrostatic and solvation properties for complex molecules

## Getting started

* Please [register your use](http://eepurl.com/by4eQr) to ensure continued support for the APBS-PDB2PQR software.
* Read our [online documentation](http://apbs-pdb2pqr.readthedocs.io/).
* Get started with the [web server](http://nbcr-222.ucsd.edu/pdb2pqr_2.1.1/).
* Download the software [following these instructions](http://apbs-pdb2pqr.readthedocs.io/en/latest/downloads.html).

-----

## A note from the developers
Dear APBS and PDB2PQR users,

We have developed a survey to better understand our user desires. In the continued development of these software, we aim to: (1) maintain and update APBS and PDB2PQR to ensure their suitability for both end-users and software developers in the biomedical community, (2) improve performance and incorporate new features based on new algorithms and user feedback, and (3) improve the robustness of the software and underlying models and reduce the likelihood of errors by inexperienced users.
 
We would truly appreciate your opinion of our software, by completing the survey, which you can access at:

<https://www.surveymonkey.com/r/APBS-PDB2PQR>
 
Thank you for your help.

Sincerely,

The APBS and PDB2PQR Development Team

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
