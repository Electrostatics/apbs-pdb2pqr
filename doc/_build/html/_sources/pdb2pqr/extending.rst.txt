Extending PDB2PQR
=================

There are several ways to extend PDB2PQR to implement new functionality and new force fields.

===================================
Custom charge and radius parameters
===================================

------------------------------------------------------------
Adding a few additional parameters to an existing forcefield
------------------------------------------------------------

If you are just adding the parameters of a few residues and atoms to an existing forcefield (e.g., AMBER), you can open the forcefield data file distributed with PDB2PQR (:file:`dat/AMBER.DAT`) directly and add your parameters.
After the parameter addition, save the force field data file with your changes.
You should also update the corresponding .names file (:file:`dat/AMBER.names`) if your added residue or atom naming scheme is different from the PDB2PQR canonical naming scheme.
See :doc:`/formats/pdb2pqr-xml-names` for more information about NAMES files.

---------------------------------
Adding an entirely new forcefield
---------------------------------

The following steps outline how to add a new force field to PDB2PQR.

You will need to generate a forcefield data file (e.g., :file:`myff.DAT`) and, if your atom naming scheme of the forcefield is different from the PDB2PQR canonical naming scheme, you will also need to provide a names files (:file:`myFF.names`).
It is recommended to build your own forcefield data and names files based on existing PDB2PQR :file:`.DAT` and :file:`.names` examples provided with PDB2PQR in the :file:`dat` directory.
See :doc:`/formats/pdb2pqr-xml-names` for more information about NAMES files.
After finishing your forcefield data file and names file, these can be used with either the command line or the web server versions of PDB2PQR.

.. todo::
   
   Provide documentation of the PDB2PQR DAT and NAMES formats.

------------------------
Adding new functionality
------------------------

PDB2PQR provides an :file:`extensions` directory that allows you to add your own code to the PDB2PQR workflow.
All functions in the extensions directory are automatically loaded into PDB2PQR as command line options using the function's name, and are called after all other steps (optimization, atom addition, parameter assignment) have been completed.
As a result any available functions are particularly useful for post-processing, or for analysis without any changes to the input structure by using the ``--clean flag``.

One of the advantages of using PDB2PQR in this fashion is the ability to use built-in PDB2PQR functions.
While a full and more detailed API can be found in the PDB2PQR ``pydoc`` documentation, some useful functions are listed below, organized by PDB2PQR module:

.. py:module:: protein

.. py:class:: Protein

   Protein objects

   .. py:method:: printAtoms(atomlist, flag)

      Print a list of atoms

   .. py:method:: getResidues()

      Return a list of residues

   .. py:method:: numResidues()

      Return the number of residues

   .. py:method:: numAtoms()

      Return the number of atoms

   .. py:method:: getAtoms()

      Return a list of atom objects

   .. py:method:: getChains()

      Return a list of chains

.. py:module:: structures

.. py:class:: Chain

   Biomolecule chain objects

   .. py:method:: getResidues()

      Return a list of residues in the chain

   .. py:method:: numResidues()

      Return the number of residues in the chain

   .. py:method:: numAtoms()

      Return the number of atoms in the chain

   .. py:method:: getAtoms()

      Return a list of atom objects in the chain

.. py:class:: Residue

   Biomolecule residue object (e.g., amino acid)

   .. py:method:: numAtoms()

      Return the number of atoms in the residue

   .. py:method:: addAtom(atom)

      Add the atom object to the residue

   .. py:method:: removeAtom(name)

      Remove a specific atom from the residue

   .. py:method:: renameAtom(old, new)

      Rename atom "old" with "new"

   .. py:method:: getAtom(name)

      Return a specific atom from the residue

   .. py:method:: hasAtom(name)

     Determine if the residue has the atom "name"

.. py:class:: Atom

   The atom of a residue

   .. py:method:: getCoords()

      Return the x/y/z coordinates of the atom

   .. py:method:: isHydrogen()

      Determine if the atom is a hydrogen or not

   .. py:method:: isBackbone()

      Determine whether the atom is from the backbone

.. py:module:: utilities

.. py:function:: getAngle(c1, c2, c3)
   
   Get the angle between the three coordinate sets

.. py:function:: getDihedral(c1, c2, c3, c4)

   Get the dihedral angle from the four coordinates

.. py:function:: distance(c1, c2)

   Return the distance between the two coordinates

.. py:function:: add(c1, c2)

   Return c1 + c2

.. py:function:: subtract(c1, c2)

   Return c1 - c2

.. py:function:: cross(c1, c2)

   Return the cross product of c1 and c2

.. py:function:: dot(c1, c2)

   Return the dot product of c1 and c2

.. py:function:: normalize(c1)

   Normalize the c1 coordinates (to unit length)

.. todo::

   Incorporate PDB2PQR Python documentation into Sphinx rather than entering it here manually.

