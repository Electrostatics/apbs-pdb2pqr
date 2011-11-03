"""
    Resinter extension

    Print interaction energy between each residue pair in the protein. 

    Author:  Kyle Monson and Emile Hogan
"""

__date__ = "21 October 2011"
__author__ = "Kyle Monson, Emile Hogan"

import extensions
from src.hydrogens import Optimize
#itertools FTW!
from itertools import product, permutations, izip, count
from collections import defaultdict
from src.hydrogens import hydrogenRoutines

_titrationSets = [['ARG','AR0'],
                  ['ASP', 'ASH'],
                  ['CYS', 'CYX'],
                  ['GLU', 'GLH'],
                  ['HSD', 'HSE', 'HSP'],
                  ['HID', 'HIE', 'HIP'],
                  ['LYN', 'LYS'],
                  ['TYR', 'TYM'],
                  ['NEUTRAL-CTERM', 'CTERM'],
                  ['NEUTRAL-NTERM', 'NTERM']]

_titrationSetsMap = defaultdict(list)

for tsSet in _titrationSets:
    for ts in tsSet:
        _titrationSetsMap[ts] = tsSet
        
#loose ends.
_titrationSetsMap['HIS'] = _titrationSetsMap['HSD']
_titrationSetsMap['CYM'] = _titrationSetsMap['CYS'] 

def addExtensionOptions(extensionGroup):
    """
        Add options.
    """
    extensionGroup.add_option('--all_combinations', dest='all_combinations', action='store_true', default=False,
                      help='Remap residues to each possible titration state and rerun resinter appending output.') 

def usage():
    return 'Print interaction energy between each residue pair in the protein to {output-path}.resinter.'

#the combinations functions allow us to iterate through all possible combinations of titration states. 
def pairwiseCombinations(initialList):
    r = [len(x) for x in initialList]
    m = min(r)
    M = max(r)
    n = len(initialList)
    
    R = set()
    
    for i in range(m):
        t = [initialList[x][i] for x in range(n)]
        R.add(tuple(t))
        
    for i in range(m,M):
        t = [initialList[x][min(r[x]-1,i)] for x in range(n)]
        R.add(tuple(t))
        
    for i in range(m):
        for j in range(n):
            for k in range(i+1,r[j]):
                prejth = [initialList[x][i] for x in range(j)]

                jth = initialList[j][k]
                
                postjth = [initialList[x][i] for x in range(j+1, n)]
                
                t = prejth + [jth] + postjth
                R.add(tuple(t))
    
    for i in range(m,M):
        for j in range(n):
            if r[i] < i:
                continue
            for k in range(i+1,r[j]):
                prejth = [initialList[x][min(r[x]-1,i)] for x in range(j)]

                jth = initialList[j][k]
                
                postjth = [initialList[x][min(r[x]-1,i)] for x in range(j+1, n)]
                
                t = prejth + [jth] + postjth
                R.add(tuple(t))
                
    return R
    
def get_residue_interaction_energy(residue1, residue2):
    """
    Returns to total energy of every atom pair between the two residues.
    
    Uses Optimize.getPairEnergy and it's donor/accepter model 
    to determine energy.
    
    residue1 - "donor" residue
    residue2 - "acceptor" residue
    
    THE RESULTS OF THIS FUNCTION ARE NOT SYMMETRIC. Swapping 
    residue1 and residue2 will not always produce the same result.
    """
    energy = 0.0
    for pair in product(residue1.getAtoms(), residue2.getAtoms()):
        energy += Optimize.getPairEnergy(pair[0], pair[1])
        
    return energy

def write_residue_interaction_energies(residues, output):
    """
    Writes out the residue interaction energy for each possible
    residue pair in the protein.
    """
    residuepairs = permutations(residues, 2)
    
    for pair in residuepairs:
        energy = get_residue_interaction_energy(pair[0], pair[1])
        output.write(str(pair[0]) + ' ' + str(pair[1]) + ' ' + str(energy) + '\n')
        
def get_residue_titration_sets(residues):
    """
    Returns all possible titration states for each residue as a list of lists.
    """
    result = []
    for residue in residues:
        result.append(_titrationSetsMap.get(residue.name, [residue.name]))
        
    return result

def process_residue_set(residueSet, routines, output, clean = False,
                                                      neutraln = False,
                                                      neutralc = False,
                                                      ligand = None,
                                                      assign_only = False,
                                                      chain = False,
                                                      debump = True,
                                                      opt = True):
    routines.write(str(residueSet)+'\n')
    
    routines.removeHydrogens()
    
    for newResidueName, oldResidue, index in izip(residueSet, routines.protein.getResidues(), count()):
        if newResidueName is None:
            continue
        
        chain = routines.protein.chainmap[oldResidue.chainID]
        chainIndex = chain.residues.index(oldResidue)
        residueAtoms = oldResidue.atoms
        
        #Create the replacement residue
        newResidue = routines.protein.createResidue(residueAtoms, newResidueName)
        
        newResidue.renameResidue(newResidueName)
        #Drop it in
        routines.protein.residues[index] = newResidue
        chain.residues[chainIndex] = newResidue
    
    #Run the meaty bits of PDB2PQR  
    routines.setTermini(neutraln, neutralc)
    routines.updateBonds()
    
    if not clean and not assign_only:
        routines.updateSSbridges()

        if debump:
            routines.debumpProtein()
            
        routines.addHydrogens()

        hydRoutines = hydrogenRoutines(routines)

        if debump:
            routines.debumpProtein()  

        if opt:
            hydRoutines.setOptimizeableHydrogens()
            hydRoutines.initializeFullOptimization()
            hydRoutines.optimizeHydrogens()
        else:
            hydRoutines.initializeWaterOptimization()
            hydRoutines.optimizeHydrogens()

        # Special for GLH/ASH, since both conformations were added
        hydRoutines.cleanup()
        
    write_residue_interaction_energies(routines.protein.getResidues(), output)        
        

def write_all_residue_interaction_energies_combinations(routines, output, options):
    """
    For every titration state combination of residue output the 
    interaction energy for all possible residue pairs. 
    """
    residueNamesList = get_residue_titration_sets(routines.protein.getResidues())
    
    routines.write("Testing the following combinations\n")
    namelist = [r.name for r in routines.protein.getResidues()]
    combinationsData = zip(namelist, residueNamesList)
    for thing in combinationsData:
        routines.write(str(thing)+'\n')
        
    count = 0
    for residueSet in pairwiseCombinations(residueNamesList):
        count += 1
        process_residue_set(residueSet, routines, output, clean = options.clean,
                                                          neutraln = options.neutraln,
                                                          neutralc = options.neutralc,
                                                          ligand = options.ligand,
                                                          assign_only = options.assign_only,
                                                          chain = options.chain,
                                                          debump = options.debump,
                                                          opt = options.opt)
    
    routines.write(str(count)+' combinations tried\n')

def create_resinter_output(routines, outfile, options, all_combinations=False):
    """
    Output the interaction energy between each possible residue pair.
    """
    routines.write("Printing residue interaction energies...\n")
    
    output = extensions.extOutputHelper(routines, outfile)
    
    if all_combinations:
        write_all_residue_interaction_energies_combinations(routines, output, options)
    else:
        write_residue_interaction_energies(routines.protein.getResidues(), output)
    

def run_extension(routines, outroot, options):
    outname = outroot + ".resinter"
    with open(outname, "w") as outfile:
        create_resinter_output(routines, outfile, options, all_combinations=options.all_combinations)