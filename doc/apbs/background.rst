Solvation model background
==========================

================
Solvation Models
================

Electrostatic and solvation models can be roughly divided into two classes ([Warshel2006]_, [Roux1999]_, [Ren2012]_) explicit solvent models that treat the solvent in atomic detail and implicit solvent models that generally replace the explicit solvent with a dielectric continuum.
Each method has its strengths and weaknesses.
While explicit solvent models offer some of the highest levels of detail, they generally require extensive sampling to converge properties of interest.
On the other hand, implicit solvent models trade detail and some accuracy for the “pre-equilibration” of solvent degrees of freedom and elimination of sampling for these degrees of freedom. Implicit solvent methods are popular for a variety of biomedical research problems.

The polar solvation energy is generally associated with a difference in charging free energies in vacuum and solvent.
A variety of implicit solvent models are available to biomedical researchers to describe polar solvation; however, the most widely-used methods are currently the Generalized Born (GB) and Poisson-Boltzmann (PB) models.
GB and related methods are very fast heuristic models for estimating the polar solvation energies of biomolecular structures and therefore are often used in high-throughput applications such as molecular dynamics simulations.
PB methods can be formally derived from more detailed theories and offer a somewhat slower, but often more accurate, method for evaluating polar solvation properties and often serve as the basis for parameterization and testing of GB methods.
Finally, unlike most GB methods, PB models provide a global solution for the electrostatic potential and field within and around a biomolecule, therefore making them uniquely suited to visualization and other structural analyses, diffusion simulations, and a number of other methods which require global electrostatic properties.

The PB equation ([Fogolari2002]_, [Lamm2003]_, [Grochowski2007]_, [Baker2005]_) is a nonlinear elliptic partial differential equation of the form shown in the figure below which is solved for the electrostatic potential.
The coefficients of this equation are directly related to the molecular structure of the system under consideration.
PB theory is approximate and, as a result, has several well-known limitations which can affect its accuracy ([Grochowski2007]_, [Netz2000]_), particularly for strongly charged systems or high salt concentrations.
However, despite these limitations, PB methods are still very important for biomolecular structural analysis, modeling, and simulation.
Furthermore, these limitations are currently being addressed through new implicit solvent models and hybrid treatments which extend the applicability of PB theory while preserving some of its computational efficiency.
There are currently examples of both types of treatments which leverage APBS ([Azuara2006]_, [Chu2007]_, [Vitalis2004]_).

.. image:: /media/pb-schematic.png

PB methods provide polar solvation energies and therefore must be complemented by non-polar solvation models to provide a complete view of biomolecular solvent-solute interactions. non-polar solvation is generally associated with the insertion of the uncharged solute into solvent. There are many non-polar solvation models available; however, work by Levy et al. [Levy2003]_ as well as our own research [Wagoner2006]_ has demonstrated the importance of non-polar implicit solvent models which include treatment of attractive solute-solvent dispersion terms.
This model has been implemented in APBS and can also be easily transformed into simpler popular non-polar models (e.g., solvent-accessible surface area).
While this model can be used separately from PB to analyze non-polar contributions to solvation energy, its preferred use is coupled to the PB equation through a geometric flow model [Chen2010]_ which treats polar and non-polar interactions in the same framework and reduces the number of user-specified empirical parameters.

==================
Further Reading
==================

.. [Azuara2006] Azuara C, Lindahl E, Koehl P, Orland H, and Delarue M, PDB_Hydro: incorporating dipolar solvents with variable density in the Poisson-Boltzmann treatment of macromolecule electrostatics. Nucleic Acids Research, 2006. 34: p. W38-W42.

.. [Baker2005] Baker NA, Biomolecular Applications of Poisson-Boltzmann Methods, in Reviews in Computational Chemistry, Lipkowitz KB, Larter R, and Cundari TR, Editors. 2005, John Wiley and Sons.

.. [Chen2010] Chen Z, Baker NA, Wei GW. Differential geometry based solvation model I: Eulerian formulation, J Comput Phys, 229, 8231-58, 2010.

.. [Chu2007] Chu VB, Bai Y, Lipfert J, Herschlag D, and Doniach S, Evaluation of Ion Binding to DNA Duplexes Using a Size-Modified Poisson-Boltzmann Theory. Biophysical Journal, 2007. 93(9): p. 3202-9.

.. [Fogolari2002] Fogolari F, Brigo A, and Molinari H, The Poisson-Boltzmann equation for biomolecular electrostatics: a tool for structural biology. Journal of Molecular Recognition, 2002. 15(6): p. 377-92.

.. [Grochowski2007] Grochowski P, lstrok A, and Trylska J, Continuum molecular electrostatics, salt effects and counterion binding. A review of the Poisson-Boltzmann theory and its modifications. Biopolymers, 2007. 89(2): p. 93-113.

.. [Lamm2003] Lamm G, The Poisson-Boltzmann Equation, in Reviews in Computational Chemistry, Lipkowitz KB, Larter R, and Cundari TR, Editors. 2003, John Wiley and Sons, Inc. p. 147-366.

.. [Levy2003] Levy RM, Zhang LY, Gallicchio E, and Felts AK, On the nonpolar hydration free energy of proteins: surface area and continuum solvent models for the solute-solvent interaction energy. Journal of the American Chemical Society, 2003. 125(31): p. 9523-30.

.. [Netz2000] Netz RR and Orland H, Beyond Poisson-Boltzmann: Fluctuation effects and correlation functions. European Physical Journal E, 2000. 1(2-3): p. 203-14.

.. [Ren2012] Ren P, Chun J, Thomas DG, Schnieders M, Marucho M, Zhang J, Baker NA, Biomolecular electrostatics and solvation: a computational perspective. Quarterly Reviews of Biophysics, 2012. 45(4): p. 427-491.

.. [Roux1999] Roux B and Simonson T, Implicit solvent models. Biophysical Chemistry, 1999. 78(1-2): p. 1-20.

.. [Vitalis2004] Vitalis A, Baker NA, McCammon JA, ISIM: A program for grand canonical Monte Carlo simulations of the ionic environment of biomolecules, Molecular Simulation, 2004, 30(1), 45-61.

.. [Wagoner2006] Wagoner JA and Baker NA, Assessing implicit models for nonpolar mean solvation forces: the importance of dispersion and volume terms. Proceedings of the National Academy of Sciences of the United States of America, 2006. 103(22): p. 8331-6.

.. [Warshel2006] Warshel A, Sharma PK, Kato M, and Parson WW, Modeling electrostatic effects in proteins. Biochimica et Biophysica Acta (BBA) - Proteins & Proteomics, 2006. 1764(11): p. 1647-76.

