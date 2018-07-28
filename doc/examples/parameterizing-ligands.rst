Parameterizing ligands
======================

This section outlines the parameterization of ligands using the PEOE_PB methods (see `Czodrowski et al, 2006 <http://dx.doi.org/10.1002/prot.21110>`_ for more information).

The PDB structure `1HPX <http://www.rcsb.org/pdb/explore.do?structureId=1hpx>`_ includes HIV-1 protease complexed with an inhibitor at 2.0 Å resolution.
HIV-1 protease has two chains; residue D25 is anionic on one chain and neutral on the other -- these titration states are important in the role of D25 as an acid in the catalytic mechanism.

-------------------
Ignoring the ligand
-------------------

If we don't want to include the ligand, then the process is straightforward:

#. From the PDB2PQR server web page, enter 1HPX into the PDB ID field.

#. Choose whichever forcefield and naming schemes you prefer.

#. Under options, be sure the "Ensure that new atoms are not rebuilt too close to existing atoms", "Optimize the hydrogen bonding network", and "Use PROPKA to assign protonation states at pH" options are selected. Choose pH 7 for your initial calculations. You can select other options as well, if interested.

#. Hit the "Submit" button.

#. Once the calculations are complete, you should see a web page with a link to the PROPKA output, a new PQR file, and warnings about the ligand KNI (since we didn't choose to parameterize it in this calculation). For comparison, you might download the the `original PDB file <http://www.pdb.org/pdb/explore.do?structureId=1HPX>`_ and compare the PDB2PQR-generated structure with the original to see where hydrogens were placed.

-------------------------
Parameterizing the ligand
-------------------------

This section outlines the parameterization of ligands using the PEOE_PB methods (see DOI:`10.1002/prot.21110 <http://dx.doi.org/10.1002/prot.21110>`_).

Ligand parameterization currently requires a :doc:`MOL2-format </formats/mol2>` representation of the ligand to provide the necessary bonding information.
MOL2-format files can be obtained through the `PRODRG web server <http://davapc1.bioch.dundee.ac.uk/cgi-bin/prodrg>`_ or some molecular modeling software packages.
PRODRG provides documentation as well as several examples on ligand preparation on its web page.

We're now ready to look at the 1HPV crystal structure from above and parameterize its ligand, KNI-272.

#. From the PDB2PQR server web page, enter ``1HPX`` into the PDB ID field.

#. Choose whichever forcefield and naming schemes you prefer.

#. Under options, be sure the "Ensure that new atoms are not rebuilt too close to existing atoms", "Optimize the hydrogen bonding network", and "Assign charges to the ligand specified in a MOL2 file" options are selected (download the :download:`ligand MOL2 file </media/LIG_1HPX.mol2>`). You can select other options as well, if interested.

#. Hit the "Submit" button.

#. Once the calculations are complete, you should see a web page with a link to the new PQR file with a warning about debumping P81 (but no warnings about ligand parameterization!). 

As a second example, we use the PDB structure `1ABF <http://www.rcsb.org/pdb/explore.do?structureId=1abf>`_ of L-arabinose binding protein in complex with a sugar ligand at 1.90 Å resolution.
To parameterize both this protein and its ligand:

#. From the PDB2PQR server web page, enter `1ABF` into the PDB ID field.

#. Choose whichever forcefield and naming schemes you prefer.

#. Under options, be sure the "Ensure that new atoms are not rebuilt too close to existing atoms", "Optimize the hydrogen bonding network", and "Assign charges to the ligand specified in a MOL2 file" options are selected (download the :download:`ligand MOL2 file </media/LIG_1ABF.mol2>`). You can select other options as well, if interested.

#. Hit the "Submit" button.

#. Once the calculations are complete, you should see a web page with a link to the new PQR file with a warning about debumping P66, K295, and K306 (but no warnings about ligand parameterization!). 
