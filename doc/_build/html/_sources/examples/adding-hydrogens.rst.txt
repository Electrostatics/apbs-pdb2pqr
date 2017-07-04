Adding hydrogens and assigning parameters
=========================================

----------------
Pick a structure
----------------

Start by choosing a PDB file to process.
Either enter the 4-character PDB ID into PDB2PQR or accession number (`1FAS <http://www.rcsb.org/pdb/explore.do?structureId=1FAS>`_ is a good starting choice) or upload your own PDB file.
Note that, if you choose to enter a 4-character PDB ID, PDB2PQR will process all recognizable chains of PDB file as it was deposited in the PDB (e.g., not the biological unit, any related transformations, etc.).

-----------------
Pick a forcefield
-----------------

For most applications, the choice is easy: PARSE.
This forcefield has been optimized for implicit solvent calculation and is probably the best choice for visualization of protein electrostatics and many common types of energetic calculations for proteins.
However, AMBER and CHARMM may be more appropriate if you are attempting to compare directly to simulations performed with those force fields, require nucleic acid support, are simulating ligands parameterized with those force fields, etc.

It is also possible to upload a user-defined forcefield (e.g., to define radii and charges for ligands or unusual residues).
Please see :doc:`/pdb2pqr/extending` for more information.

--------------------
Pick a naming scheme
--------------------

This choice is largely irrelevant to electrostatics calculations but may be important for some visualization programs.
When in doubt, choose the "Internal naming scheme" which attempts to conform to IUPAC standards.

-------------------------------------
Reconstruct missing atoms (hydrogens)
-------------------------------------

Under options, be sure the "Ensure that new atoms are not rebuilt too close to existing atoms" and "Optimize the hydrogen bonding network" options are selected.
You can select other options as well, if interested.

-----------------------------
Download and view the results
-----------------------------

Download the resulting PQR file and visualize in a molecular graphics program to examine how the hydrogens were added and how hydrogen bonds were optimized.
