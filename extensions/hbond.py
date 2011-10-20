"""
    Hbond extension

    Find all hydrogen bonds as determined by the cutoffs specified.
    Uses PDB2PQR to determine donors and acceptors, and displays
    all available bonds to file. 
    
    The original bonding parameters were an angle of 20.0, 
    distance of 3.30, and using the old method for calculating 
    distance.
    
    The original parameters for WHAT-IF output was an angle of
    90.0, distance of 3.30, and using the old method for 
    calculating distance. 

    Authors:  Todd Dolinsky, Michael J Bradley, Julie Mitchell, and Kyle Monson
"""

__date__ = "17 February 2006"
__author__ = "Todd Dolinsky, Michael J Bradley, Julie Mitchell, and Kyle Monson"
# NOTE: This extension edited and updated on 05 August 2008 by Michael J Bradley
# and again on 17 October by Kyle Monson.

# NOTE: The current defaults for hbonds used below were utilized in the following study:
#Bradley MJ, Chivers PT, Baker NA. Molecular dynamics simulation of the 
#Escherichia coli NikR protein: Equilibrium conformational fluctuations reveal 
#inter-domain allosteric communication pathways. Journal of Molecular Biology, 378, 
#1155-1173, 2008.  http://dx.doi.org/10.1016/j.jmb.2008.03.010

from src.utilities import distance, getAngle
from src.routines import Cells
from math import cos
import extensions

ANGLE_CUTOFF = 30.0       # A - D - H(D) angle
DIST_CUTOFF = 3.4         # D to A distance

def addExtensionOptions(extensionGroup):
    """
        Add options to set output type, angle cutoff, distance cutoff, and distance calculating method.
    """
    extensionGroup.add_option('--whatif', dest='whatif', action='store_true', default=False,
                              help='Change hbond output to WHAT-IF format.')
    
    extensionGroup.add_option('--angle_cutoff', dest='angle_cutoff', action='store', default=ANGLE_CUTOFF,
                              help='Angle cutoff to use when creating hbond data (default %s)' % ANGLE_CUTOFF)
    
    extensionGroup.add_option('--distance_cutoff', dest='distance_cutoff', action='store', default=DIST_CUTOFF,
                              help='Distance cutoff to use when creating hbond data (default %s)' % DIST_CUTOFF)
    
    extensionGroup.add_option('--old_distance_method', dest='old_distance_method', action='store_true', default=False,
                              help='Use distance from donor hydrogen to acceptor to calculate distance used with --distance_cutoff.')

def usage():
    return 'Print a list of hydrogen bonds to {output-path}.hbond'

#TODO: replace this with a ''.format() call.
def _residueString(residue, name):
    return '%4d %-4s (%4d  ) %s     %-4s' % \
        (residue.resSeq, residue.name, residue.resSeq, residue.chainID, name)

def create_hbond_output(routines, outfile, whatif=False, 
                                           angleCutoff=ANGLE_CUTOFF, 
                                           distanceCutoff=DIST_CUTOFF, 
                                           oldDistanceMethod=False):

    routines.write("Printing hydrogen bond list...\n")
    
    output = extensions.extOutputHelper(routines, outfile)
    
    cellsize = int(distanceCutoff + 1.0 + 1.0) 
    protein = routines.protein
    routines.setDonorsAndAcceptors()
    routines.cells = Cells(cellsize)
    routines.cells.assignCells(protein)

    for donor in protein.getAtoms():

        # Grab the list of donors
        if not donor.hdonor: 
            continue
        donorhs = []
        for bond in donor.bonds:
            if bond.isHydrogen(): 
                donorhs.append(bond)
        if donorhs == []: 
            continue

        # For each donor, grab all acceptors
        closeatoms = routines.cells.getNearCells(donor)
        for acc in closeatoms:
            if not acc.hacceptor: 
                continue
            if donor.residue == acc.residue: 
                continue
            
            #TODO: do we need to do this for plain hbond stuff?
            if whatif and (donor.residue.chainID == acc.residue.chainID): 
                continue
            
            # Do new style distance check
            if not oldDistanceMethod:
                dist = distance(donor.getCoords(), acc.getCoords())
                if dist > distanceCutoff:
                    continue
            
            for donorh in donorhs:

                # Do old style distance check       
                if oldDistanceMethod:
                    dist = distance(donorh.getCoords(), acc.getCoords())
                    if dist > distanceCutoff: 
                        continue
                    
                # Do angle check
                angle = getAngle(acc.getCoords(), donor.getCoords(), donorh.getCoords())
                if angle > angleCutoff: 
                    continue
                
                if whatif:
                    if (donor.tempFactor > 60.0): 
                        continue
                    if (acc.tempFactor > 60.0): 
                        continue
                    
                    thisBstring='B' if donor.isBackbone() else 'S'
                    thatBstring='B' if acc.isBackbone() else 'S'

                    score= (1.7/dist) * cos(angle * 3.142 / 180.0)
                    output.write(_residueString(donor.residue, donor.name))
                    output.write('-> ')
                    output.write(_residueString(acc.residue, acc.name))
                    output.write('Sym=   1 Val= %6.3lf  DA=%6.2f  DHA=%6.2f (%s-%s)\n' % 
                                 (score, dist, angle, thisBstring, thatBstring))
#                    outfile.write("%4d %-4s (%4d  ) %s     %-4s-> %4d %-4s (%4d  ) %s     %-4sSym=   1 Val= %6.3lf  DA=%6.2f  DHA=%6.2f (%s-%s)\n" % \
#                      (donor.residue.resSeq,donor.residue.name,donor.residue.resSeq, donor.residue.chainID,donor.name,acc.residue.resSeq,acc.residue.name,acc.residue.resSeq, acc.residue.chainID,acc.name, score, dist, angle, thisBstring, thatBstring)) 

                else:
                    s = "Donor: %s %s\tAcceptor: %s %s\tdist: %.2f\tAngle: %.2f\n" % \
                        (donor.residue, donor.name, acc.residue, acc.name, dist, angle)
                    output.write(s) 
                    
    routines.write("\n")

def run_extension(routines, outroot, options):
    """
        Print a list of hydrogen bonds.

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
            options:   options object 
    """
    
    outname = outroot + ".hbond"
    with open(outname, "w") as outfile:
        create_hbond_output(routines, outfile, whatif=options.whatif, 
                                               angleCutoff=options.angle_cutoff,
                                               distanceCutoff=options.distance_cutoff,
                                               oldDistanceMethod=options.old_distance_method)

