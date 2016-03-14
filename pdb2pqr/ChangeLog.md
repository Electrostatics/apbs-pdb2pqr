# These are notes for the current version of PDB2PQR
Please see http://www.poissonboltzmann.org/pdb2pqr/release-history for the complete release history

## NEW FEATURES
* Replaced the Monte Carlo method for generating titration curves with graph cut. See http://arxiv.org/abs/1507.07021

## BUG FIXES
* Added a check before calculating pKa's for large interactions energies.

## CHANGES
* The networkx library is now required for pdb2pka.

## KNOWN BUGS
* If more than one extension is run from the command line and one of the extensions modifies the protein data structure it could affect the output of the other extension. The only included extensions that exhibit this problem are resinter and newresinter.
* Running ligands and PDB2PKA at the same time is not currently supported.
* PDB2PKA currently leaks memory slowly. Small jobs will use about twice the normally required RAM (ie ~14 titratable residues will use 140MB). Big jobs will use about 5 times the normally required RAM ( 60 titratable residues will use 480MB ). We are working to fix this.
