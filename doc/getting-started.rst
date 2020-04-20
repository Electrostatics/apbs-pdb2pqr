Getting started
===============

This section gives a basic overview of APBS-PDB2PQR workflows. 
It assumes that you have `registered <http://eepurl.com/by4eQr>`_ and obtained access to the software as described in :doc:`downloads`.

The basic APBS-PDB2PQR workflow involves a few simple steps, illustrated in the figure below and enumerated as:

#. Identify your molecular structure for :ref:`pdb2pqr` by specifying a `PDB ID <http://www.pdb.org>`_ or uploading your own structure.
#. Use :ref:`pdb2pqr` to specify the titration state of the system, repair missing atoms, and assign parameters (charges and radii) to the atoms of your system.
#. Run :ref:`apbs` from within :ref:`pdb2pqr` or via the command line.
#. Visualize the results from within :ref:`pdb2pqr` or via :doc:`/other-software`.

.. image:: media/APBS-and-PDB2PQR-user-flow.png

Using PDB2PQR to prepare structures
-----------------------------------

This section outlines the basic process of adding hydrogens and assigning charge/radius parameters to an otherwise complete PDB structure.

`Fasciculin-1 <http://www.pdb.org/pdb/explore.do?structureId=1FAS>`_ is a 3-finger toxin structure available at reasonably high resolution (1.9 Å) and has all its heavy atoms present in the PDB file.
We'll use :ref:`pdb2pqr` to add hydrogens to this protein and optimize their positions.

* From the PDB2PQR server web page, enter ``1FAS`` into the PDB ID field.
* Choose whichever forcefield and naming schemes you prefer.
* Under options, be sure the "Ensure that new atoms are not rebuilt too close to existing atoms" and "Optimize the hydrogen bonding network" options are selected. You can select other options as well, if interested.
* Hit the "Submit" button.

From here, you can continue with :ref:`apbs` calculations on the server or download the file for other applications.

Note that lower-resolution structures such `1A06 <http://www.pdb.org/pdb/explore.do?structureId=1A06>`_ may cause PDB2PQR to raise errors about missing segments of the protein that cannot be reconstructed.
This behavior is expected.

