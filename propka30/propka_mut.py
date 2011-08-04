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

