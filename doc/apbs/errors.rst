Caveats and sources of error
============================

===========
Model error
===========

When performing solvation calculations using APBS, it is important to keep in mind that you are using an approximate model for solvation.
Therefore, your answers may contain errors related to approximations in the model.
Many review articles have covered the nature of these approximations, we will stress the highlights below.</p>

--------------------------
Linear dielectric response
--------------------------

The Poisson-Boltzmann equation models the solvent as a dielectric continuum that responds linearly to all applied fields.
In particular, under this model, very strong fields can induce unrealistically strong polarization in the dielectric media representing the solvent and/or the solute interior.
However, molecular solvents or solutes cannot support an infinite amount of polarization: they are limited by their density, their finite dipole moments, and their finite degree of electronic polarizability.
Therefore, the continuum model assumption of linear dielectric response can break down in situations with strong electric fields; e.g., around nucleic acids or very highly-charged proteins.

-------------------------
Local dielectric response
-------------------------

The Poisson-Boltzmann equation models the solvent as a dielectric continuum that also responds locally to all applied fields. 
n other words, under this model, the local polarization at a point x is only dependent on the field at point x.
However, molecular solvents and solutes clearly don't obey this assumption: the variety of covalent, steric, and other non-bonded intra- and inter-molecular interactions ensures that the polarization at point x is dependent on solute-field interactions in a non-vanishing neighborhood around x.
One way to limit the impact of this flawed assumption, is to model solute response as "explicitly" as possible in your continuum electrostatics problems.
In other words, rather than relying upon the continuum model to reproduce conformational relaxation or response in your solute, model such response in detail through molecular simulations or other conformational sampling.

---------------------------------------------------------
Ambiguity of dielectric interfaces and coefficient values
---------------------------------------------------------

Violation of the assumptions of linear and local dielectric response in real molecular systems leads to serious ambiguity in the definition of the dielectric coefficient in the Poisson-Boltzmann equation.
In particular, while the values for bulk solvent (i.e., far away from the solute) response are well-defined, all other values of the dielectric coefficient are ambiguous.
In general, continuum models assume a constant low-dielectric value inside the solute and the bulk solvent value outside the solute.
This assumption creates tremendous sensitivity of calculation results on the placement of the dielectric interface (usually determined by solute atomic radii) and the specific value of the internal solute dielectric.
In general, errors arising from this assumption can be minimized by using internal dielectric values that are consistent with the solute atomic radii parameterization.

--------------------------------------------------
No specific ion-solvent or ion-solute interactions
--------------------------------------------------

Most Poisson-Boltzmann models assume that ions do not interact directly with the solvent: they are charges embedded in the same dielectric material as the bulk solvent.
This assumption implies that ions experience no "desolvation" penalty as they interact with the solute surface.
Additionally, most Poisson-Boltzmann models assume that ions interaction with the solute only through electrostatic and hard-sphere steric potentials.
However, this assumption neglects some of the subtlety of ion-protein interactions; in particular, dispersive interactions that can possibly lead to some degree of ion specificity.

-----------------------
Mean field ion behavior
-----------------------

Finally, the Poisson-Boltzmann model is a "mean field" description of ionic solutions.
This means that ions only experience the average influence of other ions in the system; the model neglects fluctuations in the ionic atmosphere and correlations between the ions in solution.
Such correlations and fluctuations can be very important at high ionic charge densities; e.g., for multivalent ions, high ion concentrations, or the high-density ionic regions near highly-charged biomolecules.

====================
Parameter set errors
====================

.. todo::

   Under construction; please see https://arxiv.org/abs/1705.10035 for an initial discussion.
   Saved as issue https://github.com/Electrostatics/apbs-pdb2pqr/issues/481 

======================
Structure-based errors
======================

Electrostatics calculations can be very sensitive to errors in the structure, including:

* Misplaced atoms or sidechains

* Missing regions of biomolecular structure

* Incorrect titration state assignments

Of these errors, incorrect titration states are the most common and, often, the most problematic.
The software package PDB2PQR was created to minimize all of the above problems and we recommend its use to "pre-process" structures before electrostatics calculations.

====================
Discretization error
====================

The Poisson-Boltzmann partial differential equation must be discretized in order to be solved on a computer.
APBS discretizes the equation in spacing by evaluating the problem coefficients and solving for the electrostatic potential on a set of grid (finite difference) or mesh (finite element) points.
However, this discretization is an approximation to the actual, continuously-specified problem coefficients.
Coarser discretization of coefficients and the solution reduce the overall accuracy and introduce errors into the final potential and calculated energies.

It is very important to evaluate the sensitivity of your calculated energies to the grid spacings and lengths.
In general, it is a good idea to scan a range of grid spacings and lengths before starting a problem and choose the largest problem domain with the smallest grid spacing that gives consistent results (e.g., results that don't change as you further reduce the grid spacing).

==========================
Solver and round-off error
==========================

APBS uses iterative solvers to solve the nonlinear algebraic equations resulting from the discretized Poisson-Boltzmann equation.
Iterative solvers obtain solutions to algebraic equations which are accurate within a specified error tolerance.
Current versions of APBS use a fixed error tolerance of 10\ :sup:`-6` which implies approximately 1 part per million root-mean-squared error in calculated potentials.
Such error tolerances have been empirically observed to give good accuracy in the calculated energies obtained with APBS. 

However, it is important to note that the error in potential does not necessarily directly relate to the error in the energies calculated by APBS.
In particular, most meaningful energies are calculated as differences between energies from several calculations.
While the accuracy of each separate energy can be related to the solver error tolerance, the energy difference can only be loosely bounded by the error tolerance.

This issue is illustrated in the protein kinase ligand binding example provided with APBS as ``pka-lig`` and analyzed below.
This example demonstrates that, while the errors for each calculation remain small, the overall error in the computed energy can be very large; particularly when two different methods are compared.

.. list-table:: Sensitivity of PB energies to iterative solver error tolerance (APBS 1.2)
   :header-rows: 1

   * - Error tolerance
     - Protein energy
     - Protein energy relative error (with respect to 10\ :sup:`-12` tolerance)
     - Ligand energy
     - Ligand energy relative error (with respect to 10\ :sup:`-12` tolerance)
     - Complex energy
     - Complex energy relative error (with respect to 10\ :sup:`-12` tolerance)
     - Binding energy
     - Binding energy relative error (with respect to 10\ :sup:`-12` tolerance)
   * - 1.00E-06
     - 3.01E+05
     - 2.47E-08
     - 1.05E+04
     - 1.42E-08
     - 3.11E+05
     - 2.45E-08
     - 8.08E+00
     - 7.75E-06
   * - 1.00E-09
     - 3.01E+05
     - 3.19E-11
     - 1.05E+04
     - 1.71E-11
     - 3.11E+05
     - 2.45E-08
     - 8.08E+00
     - 2.48E-09
   * - 1.00E-12
     - 3.01E+05
     - 0.00E+00
     - 1.05E+04
     - 0.00E+00
     - 3.11E+05
     - 0.00E+00
     - 8.08E+00
     - 0.00E+00


