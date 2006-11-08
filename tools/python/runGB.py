""" Python APBS Generalized Born Implementation

    This module uses APBS to parameterize the effective Born radii of the
    Generalized Born electrostatic model and calculate the electrostatic
    solvation energy.

    Justin Xiang (jxiang@ccb.wustl.edu)
    Todd Dolinsky (todd@ccb.wustl.edu)
    Nathan Baker (baker@biochem.wustl.edu)
    Washington University in St. Louis
""" 

from apbslib import *
import sys, time, getopt
import string
import math
from sys import stdout, stderr
from math import sqrt, pow, exp, pi

__author__ = "Justin Xiang, Todd Dolinsky, Nathan Baker"
__date__ = "6 November, 2006"

Python_kb = 1.3806581e-23
Python_Na = 6.0221367e+23
Python_e0 = 8.85419e-12
Python_C = 1.602117e-19
NOSH_MAXMOL = 20
NOSH_MAXCALC = 20

class APBSError(Exception):
    """ APBSError class

        The APBSError class inherits off the Exception module and returns
        a string defining the nature of the error. 
    """
    
    def __init__(self, value):
        """
            Initialize with error message

            Parameters
                value:  Error Message (string)
        """
        self.value = value
        
    def __str__(self):
        """
            Return the error message
        """
        return `self.value`

def getHeader():
    """ Get header information about APBS
        Returns (header)
            header: Information about APBS
    """

    header = "\n\n\
    ----------------------------------------------------------------------\n\
    Adaptive Poisson-Boltzmann Solver (APBS)\n\
    Version 0.4.0\n\
    \n\
    Nathan A. Baker (baker@biochem.wustl.edu)\n\
    Dept. of Biochemistry and Molecular Biophysics\n\
    Center for Computational Biology\n\
    Washington University in St. Louis\n\
    Additional contributing authors listed in the code documentation.\n\n\
    Copyright (c) 2002-2006. Washington University in St. Louis\n\
    All Rights Reserved.\n\n\
    Portions copyright (c) 1999-2002.  University of California.\n\
    Portions copyright (c) 1995.  Michael Holst.\n\n\
    This program is free software; you can redistribute it and/or modify\n\
    it under the terms of the GNU General Public License as published by\n\
    the Free Software Foundation; either version 2 of the License, or\n\
    (at your option) any later version.\n\
    \n\
    This program is distributed in the hope that it will be useful,\n\
    but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
    GNU General Public License for more details.\n\
    \n\
    You should have received a copy of the GNU General Public License\n\
    along with this program; if not, write to the Free Software\n\
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA\n\
    \n\
    Linking APBS statically or dynamically with other modules is making a\n\
    combined work based on APBS. Thus, the terms and conditions of the GNU\n\
    General Public License cover the whole combination.\n\
    \n\
    SPECIAL GPL EXCEPTION\n\
    In addition, as a special exception, the copyright holders of APBS\n\
    give you permission to combine the APBS program with free software\n\
    programs and libraries that are released under the GNU LGPL or with\n\
    code included in releases of ISIM, Ion Simulator Interface, PMV, PyMOL\n\
    SMOL, VMD, and Vision. Such combined software may be linked with APBS and\n    redistributed together in original or modified form as mere aggregation\n\
    without requirement that the entire work be under the scope of the GNU\n\
    General Public License. This special exception permission is also extended\n    to any software listed in the SPECIAL GPL EXCEPTION clauses by the PMG,\n\
    FEtk, MC, or MALOC libraries.\n\
    \n\
    Note that people who make modified versions of APBS are not obligated\n\
    to grant this special exception for their modified versions; it is\n\
    their choice whether to do so. The GNU General Public License gives\n\
    permission to release a modified version without this exception; this\n\
    exception also makes it possible to release a modified version which\n\
    carries forward this exception.\n\
    ----------------------------------------------------------------------\n\
    \n\n"

    return header

def getUsage():
    """ Get usage information about running APBS via Python
        Returns (usage)
            usage: Text about running APBS via Python
    """
    
    usage = "\n\n\
    ----------------------------------------------------------------------\n\
    This driver program calculates electrostatic solvation energy with\n\
    Poission-Boltzmann parameterized Generailized Born model. \n\
    It is invoked as:\n\n\
      python runGB.py -i apbs.in\n\n\
    Optional arguments:\n\
      -o <output_parameter>     specifies path to output GB parameter file\n\
                                (default: GB-parameters)\n\
                                this parameter file is compatible with\n\
                                readGB.py module\n\
      -m <output_matrix>        specifies path to output GB energy matrix\n\
      -h or --help              prints this help text\n\n\
    Special input file requirements:\n\
    This program can only calculate solvation energy for 1 molecule. Also,\n\
    it is assumed that the first molecule read in from the input file is\n\
    the target molecule. Additional molecules can be included only for\n\
    other purposes, such as centering etc. The elec statments specifies\n\
    the conditions for computing self energies during the parameterization\n\
    step. \n\
    ----------------------------------------------------------------------\n\n"

    return usage

def main():
    """ Main driver for testing.  Runs APBS on given input file """
    
    # Initialize the MALOC library
    startVio()

    # Initialize variables, arrays
    com = Vcom_ctor(1)
    rank = Vcom_rank(com)
    size = Vcom_size(com)
    mgparm = MGparm()
    pbeparm = PBEparm()
    mem = Vmem_ctor("Main")
    pbe = new_pbelist(NOSH_MAXMOL)
    pmg = new_pmglist(NOSH_MAXMOL)
    pmgp = new_pmgplist(NOSH_MAXMOL)
    realCenter = double_array(3)
    totEnergy = []
    nforce = int_array(NOSH_MAXCALC)
    atomforce = new_atomforcelist(NOSH_MAXCALC)
    sdie = 0
    
    # Start the main timer
    
    main_timer_start = time.clock()

    # Check invocation
    
    stdout.write(getHeader())
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:m:", ["help"])
    except getopt.GetoptError:
        stdout.write("problem with input format\n")
        stdout.write(getUsage())
        sys.exit("Incorrect usage!\n")

    output_param = "GB-parameters.dat"
    output_matrix = ""
    for o, a in opts:
        if o in ("-h", "--help"):
            stdout.write("Printing help text...\n")
            stdout.write(getUsage())
            sys.exit()
        if o == "-i":
            input_file = a
        if o == "-o":
            output_param = a
        if o == "-m":
            output_matrix = a

    # Parse the input file
    
    nosh = NOsh_ctor(rank, size)
    try: input_file
    except:
        stdout.write("input file not initiated - check input\n")
        stdout.write(getUsage())
        sys.exit("Incorrect usage!\n")
    stdout.write("Parsing input file %s...\n" % input_file)
    if NOsh_parseInputFile(nosh, input_file) != 1:
        stderr.write("main:  Error while parsing input file.\n")
        raise APBSError, "Error while parsing input file!"

    # Load the molecules using loadMolecules routine

    alist = new_valist(NOSH_MAXMOL)
    
    if loadMolecules(nosh,alist) != 1:
        stderr.write("main:  Error while loading molecules. \n")
        raise APBSError, "Error while loading molecules!"    

    # Setup reference lists

    numAtoms = get_Valist(alist,0).number
    stdout.write("Number of atoms loaded: %d\n" % numAtoms)
    reflist = []
    radlist = []
    position = [0 for i in range(numAtoms)]
    for iatom in xrange(numAtoms):
        reflist.append(Vatom_getCharge(Valist_getAtom(get_Valist(alist, 0), iatom)))
        radlist.append(Vatom_getRadius(Valist_getAtom(get_Valist(alist, 0), iatom)))
        position[iatom] = getAtomPosition(Valist_getAtom(get_Valist(alist, 0), iatom))

    # Setup calculations

    if NOsh_setupElecCalc(nosh, alist)!=1:
        stderr.write("main: Error while setting up calculations. \n")
        raise APBSError, "Error while setting up calculations. \n"
        
    # Initialize the energy arrays
    
    for i in xrange(nosh.ncalc): totEnergy.append(0.0)
    selfEnergy = []

    # Load the dielectric maps (will be reused)
    
    dielXMap = new_gridlist(NOSH_MAXMOL)
    dielYMap = new_gridlist(NOSH_MAXMOL)
    dielZMap = new_gridlist(NOSH_MAXMOL)
    if loadDielMaps(nosh, dielXMap, dielYMap, dielZMap) != 1:
        stderr.write("Error reading dielectric maps!\n")
        raise APBSError, "Error reading dielectric maps!"
    
    kappaMap = new_gridlist(NOSH_MAXMOL)
    if loadKappaMaps(nosh, kappaMap) != 1:
        stderr.write("Error reading kappa maps!\n")
        raise APBSError, "Error reading kappa maps!"
    
    # Turn off charges

    for iatom in range(numAtoms):
        Vatom_setCharge(Valist_getAtom(get_Valist(alist, 0), iatom),0.0)

    stdout.write("-----------------------------\n")
    for iatom in range(numAtoms):

        # Turn on/off charge sequentially

        if iatom == 0:
            refCharge = reflist[iatom]
            Vatom_setCharge(Valist_getAtom(get_Valist(alist, 0), iatom),refCharge)
        else:
            refCharge = reflist[iatom]
            Vatom_setCharge(Valist_getAtom(get_Valist(alist, 0), iatom),refCharge)
            Vatom_setCharge(Valist_getAtom(get_Valist(alist, 0), iatom-1),0.0)
        
        # Load other necessary maps
    
        chargeMap = new_gridlist(NOSH_MAXMOL)
        if loadChargeMaps(nosh, chargeMap) != 1:
            stderr.write("Error reading charge maps!\n")
            raise APBSError, "Error reading charge maps!"

        # Do the calculations
	stdout.write("Parameterizaing atom %d with %d PBE calculations\n" % ((iatom+1), nosh.ncalc))
        percentage = iatom*100*numAtoms**-1
        stdout.write("Roughly %f%% complete\n" % percentage)

        energylist = [0.0 for i in range(nosh.ncalc)]
        energymatrix = []
        
        for icalc in xrange(nosh.ncalc):
            stdout.write("Details of calculation %d:\n" % (icalc+1))
            calc = NOsh_getCalc(nosh, icalc)
            mgparm = calc.mgparm
            pbeparm = calc.pbeparm

            if ((sdie == 0) & (pbeparm.sdie != 1)):
                sdie = pbeparm.sdie
            
            if calc.calctype != 0:
                stderr.write("main:  Only multigrid calculations supported!\n")
                raise APBSError, "Only multigrid calculations supported!"
            
            # Routine initMG
	
            if initMG(icalc, nosh, mgparm, pbeparm, realCenter, pbe, 
                      alist, dielXMap, dielYMap, dielZMap, kappaMap, chargeMap, 
                      pmgp, pmg) != 1:
                stderr.write("Error setting up MG calculation!\n")
                raise APBSError, "Error setting up MG calculation!"
      
            # Solve the problem : Routine solveMG
	
            thispmg = get_Vpmg(pmg,icalc)

            if solveMG(nosh, thispmg, mgparm.type) != 1:
                stderr.write("Error solving PDE! \n")
                raise APBSError, "Error Solving PDE!"

            # Set partition information : Routine setPartMG

            if setPartMG(nosh, mgparm, thispmg) != 1:
                stderr.write("Error setting partition info!\n")
                raise APBSError, "Error setting partition info!"
            
            # Get the energies - the energy for this calculation
            # (calculation number icalc) will be stored in the totEnergy array
            
            ret, totEnergy[icalc] = energyMG(nosh, icalc, thispmg, 0,
                                             totEnergy[icalc], 0.0, 0.0, 0.0)

            energies = getEnergies(thispmg, get_Valist(alist, 0))
            energylist[icalc] = 0.5*energies[iatom] # remove self interaction

            stdout.write("\n")
            
        for iprint in range(nosh.nprint):
            selfEnergy.append(returnEnergy(com, nosh, energylist, iprint))
            
        stdout.write("Self energy value: %.6E kJ/mol for atom %d\n" % (selfEnergy[-1], (iatom+1)))
        stdout.write("-----------------------------\n")
        # Clean up APBS structures for next iteration
        
        killEnergy()
        killMG(nosh, pbe, pmgp, pmg)
        killChargeMaps(nosh, chargeMap)

    # Clean up dielectirc map after parameterization
    
    killKappaMaps(nosh, kappaMap)
    killDielMaps(nosh, dielXMap, dielYMap, dielZMap)

    # SI unit conversion
    
    for i in xrange(numAtoms):
        reflist[i]=reflist[i]*Python_C
    
    # Obtain Born radii from self energies
    
    bradlist = []
    dij2 = [[0.0 for i in range(numAtoms)] for j in range(numAtoms)]
    fGB = [[0.0 for i in range(numAtoms)] for j in range(numAtoms)]
    for i in xrange(numAtoms):
        brad = -pow(reflist[i],2)*(1-1/sdie)*0.5*Python_Na/(4.0*pi*Python_e0*selfEnergy[i]*1e3)
        bradlist.append(brad)

    FILE = open(output_param, "w")
    FILE.write(str(sdie)+"\n")
    FILE.write("radii\tenergy\n")
    parameters = zip(bradlist, selfEnergy)
    for i in parameters:
        print >> FILE, "\t".join(map(str,i))
    FILE.close()
    stdout.write("GB parameters ouput in %s\n" % output_param)
    stdout.write("(1st column, Born radii; 2nd column, self energy)\n")

    for i in xrange(numAtoms):
        for j in xrange(i+1):
            
            for coord in xrange(3):
                dij2[i][j] = dij2[i][j] + pow((position[i][coord]-position[j][coord])*1e-10,2)

            d = dij2[i][j]
            bradi = bradlist[i]
            bradj = bradlist[j]
            fGB[i][j] = sqrt(d+bradi*bradj*exp(-d/(4.0*bradi*bradj)))
    
    # Calculate energy
    
    Gpol = 0.0
    for i in xrange(numAtoms):
        for j in xrange(numAtoms):
            if j <= i:
                Gpol = Gpol + reflist[i]*reflist[j]/fGB[i][j]
            else:
                Gpol = Gpol + reflist[i]*reflist[j]/fGB[j][i]
    
    Gpol = -Gpol*(1-1/sdie)*0.5*1e-3*Python_Na/(4.0*pi*Python_e0)
    
    # Print result
    
    stdout.write("\nGB Energy: %.10E kJ/mol\n" % Gpol)

    # Record data
    
    if output_matrix != "":
        FILE = open(output_matrix, "w")
        term = 0.0
        for i in range(numAtoms):
            for j in range(numAtoms):
                if j<=i:
                    term = reflist[i]*reflist[j]/fGB[i][j]
                    term = -term*(1-1/sdie)*0.5*1e-3*Python_Na/(4*pi*Python_e0)
                    FILE.write(str(term)+"\t")
                else:
                    term = reflist[i]*reflist[j]/fGB[j][i]
                    term = -term*(1-1/sdie)*0.5*1e-3*Python_Na/(4*pi*Python_e0)
                    FILE.write(str(term)+"\t")
            FILE.write("\n")
        FILE.close()
        stdout.write("Energy matrix ouput in %s\n" % output_matrix)

    # Clean up for program exit
    
    killMolecules(nosh, alist)
    
    delete_Nosh(nosh)
        
    # Clean up Python structures
        
    delete_double_array(realCenter)
    delete_int_array(nforce)
    delete_atomforcelist(atomforce)
    delete_valist(alist)
    delete_gridlist(dielXMap)
    delete_gridlist(dielYMap)
    delete_gridlist(dielZMap)
    delete_gridlist(kappaMap)
    delete_gridlist(chargeMap)
    delete_pmglist(pmg)
    delete_pmgplist(pmgp)
    delete_pbelist(pbe)
    
    
    # Clean up MALOC structures
    delete_Com(com)
    delete_Mem(mem)
    stdout.write("\n")
    stdout.write("Thanks for using APBS!\n\n")

    # Stop the main timer
    main_timer_stop = time.clock()
    stdout.write("Total execution time:  %1.6e sec\n" % (main_timer_stop - main_timer_start))

 
if __name__ == "__main__": main()
