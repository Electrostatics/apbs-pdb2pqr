"""
    Resinter extension

    Print interaction energy between each residue pair in the protein. 

    Author:  Kyle Monson and Emile Hogan
"""

__date__ = "21 October 2011"
__author__ = "Kyle Monson and Emile Hogan"

import extensions
from src.hydrogens import Optimize
#itertools FTW!
from itertools import product, permutations

def usage():
    return 'Print interaction energy between each residue pair in the protein to {output-path}.resinter.'

def get_residue_interaction_energy(residue1, residue2):
    energy = 0.0
    for pair in product(residue1.getAtoms(), residue2.getAtoms()):
        energy += Optimize.getPairEnergy(pair[0], pair[1])
        
    return energy

def create_resinter_output(routines, outfile):
    routines.write("Printing residue interaction energies...\n")
    
    output = extensions.extOutputHelper(routines, outfile)
    
    residuepairs = permutations(routines.protein.getResidues(), 2)
    
    for pair in residuepairs:
        energy = get_residue_interaction_energy(pair[0], pair[1])
        output.write(str(pair[0]) + ' ' + str(pair[1]) + ' ' + str(energy) + '\n')

def run_extension(routines, outroot, options):
    outname = outroot + ".resinter"
    with open(outname, "w") as outfile:
        create_resinter_output(routines, outfile)