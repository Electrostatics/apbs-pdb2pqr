from apbslib import *
import sys, time

from sys import stdout, stderr

class APBSError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return `self.value`

Python_kb = 1.3806581e-23
Python_Na = 6.0221367e+23
NOSH_MAXMOL = 20
NOSH_MAXCALC = 20

header = "\n\n\
    ----------------------------------------------------------------------\n\
    Adaptive Poisson-Boltzmann Solver (APBS)\n\
    Version 0.2.5 (November 21, 2002)\n\
    \n\
    Nathan A. Baker (baker@biochem.wustl.edu)\n\
    Dept. of Biochemistry and Molecular Biophysics\n\
    Center for Computational Biology\n\
    Washington University in St. Louis\n\
    Additional contributing authors listed in the code documentation.\n\n\
    Copyright (c) 2003. Washington University in St. Louis\n\
    All Rights Reserved.\n\n\
    Portions copyright (c) 1999-2002.  University of California.\n\
    Portions copyright (c) 1995.  Michael Holst.\n\n\
    Permission to use, copy, modify, and distribute this software and its\n\
    documentation for educational, research, and not-for-profit purposes,\n\
    without fee and without a signed licensing agreement, is hereby granted,\n\
    provided that the above copyright notice, this paragraph and the\n\
    following two paragraphs appear in all copies, modifications, and\n\
    distributions.\n\n\
    IN NO EVENT SHALL THE AUTHORS BE LIABLE TO ANY PARTY FOR DIRECT,\n\
    INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST\n\
    PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,\n\
    EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH\n\
    DAMAGE.\n\n\
    THE AUTHORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT\n\
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n\
    PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF\n\
    ANY, PROVIDED HEREUNDER IS PROVIDED \"AS IS\".  THE AUTHORS HAVE NO\n\
    OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR\n\
    MODIFICATIONS.\n\
    ----------------------------------------------------------------------\n\
    \n\n"

usage = "\n\n\
    ----------------------------------------------------------------------\n\
    This driver program calculates electrostatic potentials, energies,\n\
    and forces using both multigrid and finite element methods.\n\
    It is invoked as:\n\n\
      apbs apbs.in\n\n\
    where apbs.in is a formatted input file.\n\
    ----------------------------------------------------------------------\n\n"

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
totEnergy = double_array(NOSH_MAXCALC)
qfEnergy = double_array(NOSH_MAXCALC)
qmEnergy = double_array(NOSH_MAXCALC)
dielEnergy = double_array(NOSH_MAXCALC)
npEnergy = double_array(NOSH_MAXCALC)
nenergy = int_array(NOSH_MAXCALC)
nforce = int_array(NOSH_MAXCALC)
atomforce = new_atomforcelist(NOSH_MAXCALC)

# Start the main timer
main_timer_start = time.clock()

# Check invocation
stdout.write(header)
if len(sys.argv) != 2:
	stderr.write("main:  Called with %d arguments!\n" % len(sys.argv))
	stderr.write(usage)
	raise APBSError, "Oops!"

# Parse the input file
nosh = NOsh()
NOsh_ctor2(nosh, rank, size)
input_file = sys.argv[1]
stdout.write("Parsing input file %s...\n" % input_file)
if NOsh_parseFile(nosh, input_file) != 1:
	stderr.write("main:  Error while parsing input file.\n")
	raise APBSError, "Oops!"

# Load the molecules using loadMolecules routine

alist = new_valist(NOSH_MAXMOL)
if loadMolecules(nosh,alist) != 1:
	  stderr.write("main:  Error while loading molecules. \n")
	  raise APBSError, "Oops!"

# Load the dieletric maps

dielXMap = new_gridlist(NOSH_MAXMOL)
dielYMap = new_gridlist(NOSH_MAXMOL)
dielZMap = new_gridlist(NOSH_MAXMOL)
if loadDielMaps(nosh, dielXMap, dielYMap, dielZMap) != 1:
	stderr.write("Error reading dielectric maps!\n")
	raise APBSError, "Oops!"

# Load the kappa maps
kappaMap = new_gridlist(NOSH_MAXMOL)
if loadKappaMaps(nosh, kappaMap) != 1:
	stderr.write("Error reading kappa maps!\n")
	raise APBSError, "Oops!"

# Load the charge maps
chargeMap = new_gridlist(NOSH_MAXMOL)
if loadChargeMaps(nosh, chargeMap) != 1:
	stderr.write("Error reading charge maps!\n")
	raise APBSError, "Oops!"

# Do the calculations

stdout.write("Preparing to run %d PBE calculations. \n" % nosh.ncalc)

for icalc in xrange(nosh.ncalc):
	stdout.write("---------------------------------------------\n")
	calc = NOsh_getCalc(nosh, icalc)
	mgparm = calc.mgparm
	pbeparm = calc.pbeparm
	if calc.calctype != 0:
		stderr.write("main:  Only multigrid calculations supported!\n")
		raise APBSError, "Oops!"
	stdout.write("CALCULATION #%d:  MULTIGRID\n" % (icalc+1))
	stdout.write("Setting up problem...\n")
	
	# Routine initMG
	
	if initMG(icalc, nosh, mgparm, pbeparm, realCenter, pbe, 
              alist, dielXMap, dielYMap, dielZMap, kappaMap, chargeMap, 
              pmgp, pmg) != 1:
		stderr.write("Error setting up MG calculation!\n")
		raise APBSError, "Oops!"
	
	# Print problem parameters 
	
	printMGPARM(mgparm, realCenter)
	printPBEPARM(pbeparm)
	Python_kbT = Python_kb*Python_Na*(pbeparm.temp)/1000
	
	# Solve the problem : Routine solveMG
	
	thispmg = get_Vpmg(pmg,icalc)

	if solveMG(nosh, thispmg, mgparm.type) != 1:
		stderr.write("Error solving PDE! \n")
		raise APBSError, "Oops!"

	# Set partition information : Routine setPartMG

	if setPartMG(nosh, mgparm, thispmg) != 1:
		stderr.write("Error setting partition info!\n")
		raise APBSError, "Oops"
	
	# Write out energies : Routine energyMG, npenergyMG
	# Create pointers to variables that store energies,
	#   place in appropriate array, and delete pointers
	# See SWIG documentation for pointer function info
	
	totEng = ptrcreate("double",0.0)
	qfEng = ptrcreate("double",0.0)
	qmEng = ptrcreate("double",0.0)
	dielEng = ptrcreate("double",0.0)
	npEng = ptrcreate("double",0.0)
	neng = ptrcreate("int",0.0)

	energyMG(nosh, icalc, thispmg, nenergy,  totEng, qfEng, qmEng, dielEng)
	#npenergyMG(nosh, icalc, thispmg, nenergy, npEng)
	
	ptrset(totEnergy,ptrvalue(totEng),icalc)
	ptrset(qfEnergy,ptrvalue(qfEng),icalc)
	ptrset(qmEnergy,ptrvalue(qmEng),icalc)
	ptrset(dielEnergy,ptrvalue(dielEng),icalc)
	ptrset(npEnergy,ptrvalue(npEng),icalc)
	ptrset(nenergy,ptrvalue(neng),icalc)
	ptrfree(totEng)
	ptrfree(qfEng)
	ptrfree(qmEng)
	ptrfree(dielEng)
	ptrfree(npEng)
	ptrfree(neng)
	
	# Set partition information

	nfor = ptrcreate("int",0.0)
	forceMG(mem, nosh, pbeparm, thispmg, nfor, atomforce, alist)
	ptrset(nforce,ptrvalue(nfor), icalc)
	ptrfree(nfor)
	
	# Write out data from MG calculations : Routine writedataMG	
	writedataMG(rank, nosh, pbeparm, thispmg)
	
	# Write out matrix from MG calculations	
	writematMG(rank, nosh, pbeparm, thispmg)

# Handle print statements

if nosh.nprint > 0:
	stdout.write("---------------------------------------------\n")
	stdout.write("PRINT STATEMENTS\n")
for iprint in xrange(nosh.nprint):
	if NOsh_printWhat(nosh, iprint) == 0:
		printEnergy(com, nosh, totEnergy, iprint)
	elif NOsh_printWhat(nosh, iprint) == 1:
		printForce(com, nosh, nforce, atomforce, iprint)
	else:
		stdout.write("Undefined PRINT keyword!\n")
		break
	
stdout.write("----------------------------------------\n")
stdout.write("CLEANING UP AND SHUTTING DOWN...\n")

# Clean up APBS structures
killForce(mem, nosh, nforce, atomforce)
killEnergy()
killMG(nosh, pbe, pmgp, pmg)
killChargeMaps(nosh, chargeMap)
killKappaMaps(nosh, kappaMap)
killDielMaps(nosh, dielXMap, dielYMap, dielZMap)
killMolecules(nosh, alist)
#NOsh_dtor(getNosh(nosh))

# Clean up MALOC structures
Vcom_dtor(getCom(com))
Vmem_dtor(getMem(mem))
stdout.write("\n")
stdout.write("Thanks for using APBS!\n\n")

# Stop the main timer
main_timer_stop = time.clock()
stdout.write("Total execution time:  %1.6e sec\n" % 
  (main_timer_stop - main_timer_start))

