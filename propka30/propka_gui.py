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
import Source.version as propka
import Source.lib as lib
from Source.protein import Protein
from Source.mutate import makeCompositeAtomsDictionary
 


def main():
    """
    This is a tailor-made propka-executable for propka's GUI, but also useful for commandline execution
    """

    # I. preliminaries

    # loading options, flaggs and arguments; making a version object
    options, pdbfiles = lib.loadOptions()
    version = propka.makeVersion(label=options.version_label)

    # creating the protein object
    myProtein = Protein(pdbfile=pdbfiles[0], options=options)

    # creating a dictionary with atom objects needed for e.g. alignment mutations
    atoms = makeCompositeAtomsDictionary(protein=myProtein, pdbfiles=pdbfiles, options=options)


    # II. optimise 'single-site' mutation and selecting the good determinants
    mutations = []
    for mutation in options.mutations:
        # II.a. combinatorial search of 'determinants' for each site
        best_mutation = myProtein.optimizeMutationDeterminants(atoms=atoms, mutation=mutation, version=version, options=options)
        if best_mutation != None:
            mutations.append( best_mutation )

    # III. combinatorial search of single-site mutations with their resulting 'determinants'
    if len(mutations) > 0:
      best_mutation = myProtein.optimizeMultipleMutations(mutations=mutations, atoms=atoms, version=version, options=options)
    else:
      print("Could not find any mutation combination more stable than WT\n")
      sys.exit(8)



if __name__ == '__main__': main()

