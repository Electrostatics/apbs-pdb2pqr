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
from Source.protein import Protein
 


def main():
    """
    Simple main that creates the proteins, calculates pKa values, and prints pKa files
    """

    # loading options, flaggs and arguments
    options, pdbfiles = lib.loadOptions()

    for pdbfile in pdbfiles:

      # creating protein object
      myProtein = Protein(pdbfile=pdbfile, options=options)

      # calculating pKa values for ionizable residues
      myProtein.calculatePKA(options=options)
      # printing pka file
      myProtein.writePKA(options=options)


if __name__ == '__main__': main()

