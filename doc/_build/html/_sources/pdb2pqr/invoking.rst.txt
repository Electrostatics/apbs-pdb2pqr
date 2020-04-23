Invoking PDB2PQR
================

.. note::

   Please `register <http://eepurl.com/by4eQr>`_ before using PDB2PQR!

Most users will use PDB2PQR through `the web server <http://apbs.poissonboltzmann.org/>`_ (after `registering <http://eepurl.com/by4eQr>`_, of course).
However, it is also possible to install local versions of PDB2PQR and run these through the command line.
This version of the software offers an expanded range of options and can also be customized with user extensions.

.. todo::

   Make sure user extensions are adequately documented.

The command line PDB2PQR is invoked as

.. code-block:: bash

   $ python pdb2pqr.py [options] --ff={forcefield} {path} {output-path}

This module takes a PDB file as input and performs optimizations before yielding a new PQR-style file in ``{output-path}``.
If ``{path}`` is a `PDB ID <http://www.rcsb.org/pdb/staticHelp.do?p=help/advancedsearch/pdbIDs.html>`_ it will automatically be retrieved from the online PDB archive.

=======
Options
=======

.. todo::

   Replace this section with automagically generated documentation ala https://sphinx-argparse.readthedocs.io/en/stable/

-----------------
Mandatory options
-----------------

One of the options must be used for all PDB2PQR runs:

``--ff=FIELD_NAME``
  The forcefield to use - currently amber, charmm, parse, tyl06, peoepb and swanson are supported.

``--userff=USER_FIELD_FILE``
  A user-created forcefield file. Requires ``--usernames`` and overrides ``--ff``.

``--clean``
  Do no optimization, atom addition, or parameter assignment, just return the original PDB file in aligned format.
  Overrides ``--ff`` and ``--userff`` options.

-----------------------------------
Titration state calculation options
-----------------------------------

``--ph-calc-method=PH_METHOD``
  Method used to calculate ph values.
  If a pH calculation method is selected, pKa values will be calculated and titratable residues potentially modified after comparison with the pH value supplied by ``--with_ph`` for each titratable residue

  ``PROPKA``
    Use PROPKA to calculate pKa values.
    PROPKA options include:

    ``--propka-reference=PROPKA_REFERENCE``
      Setting which reference values to use for stability calculations.
      See PROPKA 3.0 documentation for details.

    ``--propka-verbose``
      Print extra-verbose PROPKA information to :file:`stdout`.

  ``PDB2PKA``
    Use PDB2PKA to calculate pKa values.
    Requires the use of the PARSE force field.
    Large proteins can take a very long time to run using this method.
    PDB2PKA options include:

    ``--pdb2pka-out=PDB2PKA_OUT``
      Output directory for PDB2PKA results. Defaults to :file:`pdb2pka_output`.

    ``--pdb2pka-resume``
      Resume run from state saved in output directory.

    ``--pdie=PDB2PKA_PDIE``
      Protein dielectric constant.
      Defaults to 8.

    ``--sdie=PDB2PKA_SDIE``
      Solvent dielectric constant.
      Defaults to 80.

    ``--pairene=PDB2PKA_PAIRENE``
      Cutoff energy in (kT) for calculating interaction energies.
      Default: 1.0.

``--with-ph=PH``
  pH values to use when applying the results of the selected pKa calculation method to assign titration states.
  Defaults to 7.0.

.. todo::

   Make sure that PDB2PQR force fields are referenced and described.

---------------
General options
---------------

``--apbs-input``
  Create a template APBS input file based on the generated PQR file.
  Also create a Python pickle for using these parameters in other programs.

``--assign-only``
  Only assign charges and radii - do not add atoms, debump, or optimize.

``--chain``
  Keep the PDB chain ID in the output PQR file.

``--drop-water``
  Drop waters before processing protein.
  Currently recognized and deleted are the following water types: ``HOH``, ``WAT``

``--ffout=FIELD_NAME``
  Instead of using the standard canonical naming scheme for residue and atom names, use the names from the given forcefield.
  Currently amber, charmm, parse, tyl06, peoepb and swanson are supported.

``-h`` or ``--help``
  Print help message and exit.

``--ligand=PATH``
  Calculate the parameters for the ligand in mol2 format at the given path.
  PDB2PKA must be compiled.

``--neutraln``
  Make the N-terminus of this protein neutral (default charge state is determined by pH and pKa).
  Requires PARSE force field.

``--neutralc``
  Make the C-terminus of this protein neutral (default is charged).
  Requires PARSE force field.

``--nodebump``
  Do not perform the debumping operation to remove steric clashes.
  See debumping_ for more information.

``--noopts``
  Do not perform hydrogen bond optimization.
  See hbondopt_ for more information.

``--typemap``
  Create a map of atom types in the molecule.

``--usernames=USER_NAME_FILE``
  The user created names file to use. Required if using ``--userff``.

``-v`` or ``--verbose``
  Print additional information to stdout.

``--version``
  Show program's version number and exit

``--whitespace``
  Insert whitespaces between atom name and residue name, between x and y, and between y and z.

``--include_header``
  Include pdb header in pqr file.

  .. warning::

     The resulting PQR file will not with APBS versions prior to 1.5.

-----------------
Extension options
-----------------

These options refer to extension modules distributed with PDB2PQR.

``--chi``
  Print the per-residue sidechain chi angles to :file:`{output-path}.chi`

``--contact``
  Print a list of contacts to :file:`{output-path}.con`

``--hbond``
  Print a list of hydrogen bonds to :file:`{output-path}.hbond`.
  Additional options for this extension include:

  ``--whatif``
    Change hbond output to WHAT-IF format.

  ``--angle_cutoff=ANGLE_CUTOFF``
    Angle cutoff to use when creating hbond data (default 30.0 deg)

  ``--distance_cutoff=DISTANCE_CUTOFF``
    Distance cutoff to use when creating hbond data (default 0.34 nm)

  ``--old_distance_method``
    Use distance from donor hydrogen to acceptor to calculate distance used with ``--distance_cutoff``.

``--rama``
  Print the per-residue phi and psi angles to :file:`{output-path}.rama` for Ramachandran plots.
  Options for this extension include:

  ``--phi_only``
    Only include phi angles in output. Rename output file :file:`{output-path}.phi`

  ``--psi_only``
    Only include psi angles in output. Rename output file :file:`{output-path}.psi`

``--resinter`` or ``--newresinter``
  Print interaction energy between each residue pair in the protein to :file:`{output-path}.resinter` or :file:`{output-path}.newresinter`
  Additional options for this extension include:

  ``--residue_combinations``
    Remap residues to different titration states and rerun ``--resinter``, appending the output.
    Consider only the minimum number of whole protein titration combinations needed to test each possible pairing of residue titration states.
    Normally used with ``--noopt``.
    If a protein titration state combination results in a pair of residue being re-tested in the same individual titration states, a warning will be generated if the re-tested result is different.
    This warning should not be possible if used with ``--noopt``.

  ``--all_residue_combinations``
    Remap residues to ALL possible titration state combinations and rerun resinter appending output.
    Results with ``--noopt`` should be the same as ``--residue_combinations``.
    Runs considerably slower than ``--residue_combinations`` and generates the same type of warnings.
    Use without ``--noopt`` to discover how hydrogen optimization affects residue interaction energies via the warnings in the output.

``--salt``
  Print a list of salt bridges to :file:`{output-path}.salt`

``--summary``
  Print protein summary information to :file:`{output-path}.summary`

===================
Method descriptions
===================

.. _debumping:

---------
Debumping
---------

Unless otherwise instructed with ``--nodebump``, PDB2PQR will attempt to remove steric clashes (debump) between residues.

To determine if a residue needs to be debumped, PDB2PQR compares its atoms to all nearby atoms.
With the exception of donor/acceptor pairs and CYS residue SS bonded pairs, a residue needs to be debumped if any of its atoms are within cutoff distance of any other atoms.
The cutoff is 1.0 angstrom for hydrogen/hydrogen collisions, 1.5 angstrom for hydrogen/heavy collisions, and 2.0 angstrom otherwise. 

Considering the atoms that are conflicted, PDB2PQR changes selected dihedral angle configurations in increments of 5.0 degrees, looking for positions where the residue does not conflict with other atoms.
If modifying a dihedral angle does not result in a debumped configuration then the dihedral angle is reset and the next one is tried.
If 10 angles are tried without success the algorithm reports failure. 

.. warning::

   It should be noted that this is not an optimal solution.
   This method is not guaranteed to find a solution if it exists and will accept the first completely debumped state found, not the optimal state. 

   Additionally, PDB2PQR does not consider water atoms when looking for conflicts.

.. _hbondopt:

--------------------------
Hydrogen bond optimization
--------------------------

Unless otherwise indicated with ``--noopts``, PDB2PQR will attempt to add hydrogens in a way that optimizes hydrogen bonding.

The hydrogen bonding network optimization seeks, as the name suggests, to optimize the hydrogen bonding network of the protein.
Currently this entails manipulating the following residues:

* Flipping the side chains of HIS (including user defined HIS states), ASN, and GLN residues;
* Rotating the sidechain hydrogen on SER, THR, TYR, and CYS (if available);
* Determining the best placement for the sidechain hydrogen on neutral HIS, protonated GLU, and protonated ASP;
* Optimizing all water hydrogens.

------------------------------------
Titration states
------------------------------------

PDB2PQR has the ability to recognize certain protonation states and keep them fixed during optimization.
To use this feature manually rename the residue name in the PDB file as follows:

Neutral ASP
  ASH

Negative CYS:
  CYM

Neutral GLU:
  GLH

Neutral HIS:
  HIE or HSE (epsilon-protonated); HID or HSD (delta-protonated)

Positive HIS:
  HIP or HSP

Neutral LYS
  LYN

Negative TYR
  TYM

PDB2PQR is unable to assign charges and radii when they are not available in the forcefield - thus this warning message will occur for most ligands unless a MOL2 file is provided for the ligand with the ``--ligand`` option.
Occasionally this message will occur in error for a standard amino acid residue where an atom or residue may be misnamed.
However, some of the protonation states derived from the PROPKA results are not supported in the requested forcefield and thus PDB2PQR is unable to get charges and radii for that state.
PDB2PQR currently supports the following states as derived from PROPKA:

================== ============= ============== =============
Protonation State  AMBER Support CHARMM Support PARSE Support
================== ============= ============== =============
Neutral N-Terminus No            No             Yes
Neutral C-Terminus No            No             Yes
Neutral ARG        No            No             No
Neutral ASP        Yes [#but]_   Yes            Yes
Negative CYS       Yes [#but]_   No             Yes
Neutral GLU        Yes [#but]_   Yes            Yes
Neutral HIS        Yes           Yes            Yes
Neutral LYS        Yes [#but]_   No             Yes
Negative TYR       No            No             Yes
================== ============= ============== =============

.. [#but] Only if residue is not a terminal residue; if the residue is terminal it will not be set to this state.
