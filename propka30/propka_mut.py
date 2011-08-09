#!/usr/bin/env python
#
# * This library is free software; you can redistribute it and/or
# * modify it under the terms of the GNU Lesser General Public
# * License as published by the Free Software Foundation; either
# * version 2.1 of the License, or (at your option) any later version.
# *
# * This library is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# * Lesser General Public License for more details.
#

#propka3.0, revision 182                                                                      2011-08-09
#-------------------------------------------------------------------------------------------------------
#--                                                                                                   --
#--                                   PROPKA: A PROTEIN PKA PREDICTOR                                 --
#--                                                                                                   --
#--                                VERSION 1.0,  04/25/2004, IOWA CITY                                --
#--                                             BY HUI LI                                             --
#--                                                                                                   --
#--                               VERSION 2.0,  11/05/2007, IOWA CITY/COPENHAGEN                      --
#--                                BY DELPHINE C. BAS AND DAVID M. ROGERS                             --
#--                                                                                                   --
#--                              VERSION 3.0,  xx/xx/2010, COPENHAGEN                                 --
#--                              BY MATS H.M. OLSSON AND CHRESTEN R. SONDERGARD                       --
#--                                                                                                   --
#-------------------------------------------------------------------------------------------------------
#
#
#-------------------------------------------------------------------------------------------------------
# References:
#
#   Very Fast Empirical Prediction and Rationalization of Protein pKa Values
#   Hui Li, Andrew D. Robertson and Jan H. Jensen
#   PROTEINS: Structure, Function, and Bioinformatics 61:704-721 (2005)
#
#   Very Fast Prediction and Rationalization of pKa Values for Protein-Ligand Complexes
#   Delphine C. Bas, David M. Rogers and Jan H. Jensen
#   PROTEINS: Structure, Function, and Bioinformatics 73:765-783 (2008)
#
#   PROPKA3: Consistent Treatment of Internal and Surface Residues in Empirical pKa predictions
#   Mats H.M. Olsson, Chresten R. Sondergard, Michal Rostkowski, and Jan H. Jensen
#   Journal of Chemical Theory and Computation, to be submitted (2010)
#-------------------------------------------------------------------------------------------------------

import string, re, sys, os, math

import Source.lib as lib
import Source.mutate as mutate
from Source.protein import Protein
from Source.pdb import readPDB
import Source.version as propka
 


def main():
    """
    Simple main that creates the proteins, mutates the protein, calculates pKa values, and prints pKa files
    """

    # I. preliminaries

    # loading options, flaggs and arguments
    options, pdbfiles = lib.loadOptions()
    version = propka.makeVersion(label=options.version_label)

    # creating protein object and calculating reference folding energy
    myProtein = Protein(pdbfile=pdbfiles[0], options=options)
    myProtein.calculatePKA(version=version, options=options)
    dG_ref = myProtein.calculateFoldingEnergy(options=options)

    # creating a dictionary with atom objects needed for e.g. alignment mutations
    atoms = mutate.makeCompositeAtomsDictionary(protein=myProtein, pdbfiles=pdbfiles, options=options)


    # II. making mutations and calculating the folding energy
    for mutation in options.mutations:
      print(mutation)
      newProtein = mutate.makeMutatedProtein(myProtein, atoms=atoms, mutation=mutation, options=options)

      # calculating pKa values for ionizable residues
      newProtein.calculatePKA(version=version, options=options)
      newProtein.writePKA(options=options)
      dG_mut = newProtein.calculateFoldingEnergy(options=options)
      print("staliblization:")
      print("%8.2lf %6.2lf  %s" % (dG_ref, dG_ref-dG_ref, "WT"))
      print("%8.2lf %6.2lf  %s" % (dG_mut, dG_mut-dG_ref, "Mutant"))
      newProtein.writePDB(hydrogens=True, options=options)


if __name__ == '__main__': main()

