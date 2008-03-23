""" Python APBS Generalized Born Implementation

    This module uses APBS to parameterize the effective Born radii of the
    Generalized Born electrostatic model and calculate the electrostatic
    solvation energy.

    Justin Xiang (jxiang@ccb.wustl.edu)
    Todd Dolinsky (todd@ccb.wustl.edu)
    Nathan Baker (baker@biochem.wustl.edu)
    Washington University in St. Louis

	APBS -- Adaptive Poisson-Boltzmann Solver

	  Nathan A. Baker (baker@biochem.wustl.edu)
	  Dept. Biochemistry and Molecular Biophysics
	  Center for Computational Biology
	  Washington University in St. Louis

	  Additional contributing authors listed in the code documentation.

	Copyright (c) 2002-2007, Washington University in St. Louis.
	Portions Copyright (c) 2002-2007.  Nathan A. Baker
	Portions Copyright (c) 1999-2002.  The Regents of the University of California.
	Portions Copyright (c) 1995.  Michael Holst

	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met: 

	* Redistributions of source code must retain the above copyright notice, this
	list of conditions and the following disclaimer.  

	* Redistributions in binary form must reproduce the above copyright notice,
	this list of conditions and the following disclaimer in the documentation
	and/or other materials provided with the distribution.

	* Neither the name of Washington University in St. Louis nor the names of its
	contributors may be used to endorse or promote products derived from this
	software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
	A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
	EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
	PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
	PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
""" 

import sys, time, getopt
import string
import math
from sys import stdout, stderr
from math import sqrt, pow, exp, pi

__author__ = "Todd Dolinsky, Nathan Baker"
__date__ = "July 2007"

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
    Version 0.5.1\n\
    \n\
    APBS -- Adaptive Poisson-Boltzmann Solver\n\
    \n\
    Nathan A. Baker (baker@biochem.wustl.edu)\n\
    Dept. Biochemistry and Molecular Biophysics\n\
    Center for Computational Biology\n\
    Washington University in St. Louis\n\
    \n\
    Additional contributing authors listed in the code documentation.\n\
    \n\
    Copyright (c) 2002-2007, Washington University in St. Louis.\n\
    Portions Copyright (c) 2002-2007.  Nathan A. Baker\n\
    Portions Copyright (c) 1999-2002.  The Regents of the University of California.\n\
    Portions Copyright (c) 1995.  Michael Holst\n\
    \n\
    All rights reserved.\n\
    \n\
    Redistribution and use in source and binary forms, with or without\n\
    modification, are permitted provided that the following conditions are met:\n\
    \n\
    * Redistributions of source code must retain the above copyright notice, this\n\
      list of conditions and the following disclaimer.\n\
    \n\
    * Redistributions in binary form must reproduce the above copyright notice,\n\
      this list of conditions and the following disclaimer in the documentation\n\
      and/or other materials provided with the distribution.\n\
    \n\
    * Neither the name of Washington University in St. Louis nor the names of its\n\
      contributors may be used to endorse or promote products derived from this\n\
      software without specific prior written permission.\n\
    \n\
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n\
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n\
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n\
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR\n\
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,\n\
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,\n\
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR\n\
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF\n\
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING\n\
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\n\
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\
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
    This driver program calculates electrostatic solvation energy from the\n\
    Generalized Born model based on a provided parameter file.\n\
    It is invoked as:\n\n\
      python readGB.py -p parameter_file -i structure.pqr\n\n\
    The parameter_file can be generated using runGB modules. It should\n\
    contain either, or both of the column fields: Effective Born radii and\n\
    self energy. The first line of the parameter_file should indicate the\n\
    value of solvent dielectric. The second line specifies what the column\n\
    fields are, by flagging 'radii' and/or 'energy'.\n\n\
    Optional arguments:\n\
      -o <output_parameter>     specifies path to output GB parameter file\n\
      -m <output_matrix>        specifies path to output GB energy matrix\n\
      -h or --help              prints this help text\n\
    ----------------------------------------------------------------------\n\n"

    return usage

def main():

    main_timer_start = time.clock()
    
    sdie = 0.0

    # Method: 1 = from energy; 0 = from radii
    method = 0
    column = -1
    
    # Check invocation
    stdout.write(getHeader())
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hp:i:o:m:", ["help"])
    except getopt.GetoptError:
        stdout.write("problem with input format\n")
        stdout.write(getUsage())
        sys.exit("Incorrect usage!\n")
        
    output_param = ""
    output_matrix = ""
    for o, a in opts:
        if o in ("-h", "--help"):
            stdout.write("Printing help text...\n")
            stdout.write(getUsage())
            sys.exit()
        if o == "-p":
            if a == "":
                stdout.write("parameter file not specified\n")
                stdout.write(getUsage())
                sys.exit("Incorrect usage!\n")
            else:
                param_file = a
        if o == "-i":
            if a == "":
                stdout.write("input pqr file not specified\n")
                stdout.write(getUsage())
                sys.exit("Incorrect usage!\n")
            else:
                pqr_file = a
        if o == "-o":
            output_param = a
        if o == "-m":
            output_matrix = a

    # Parse the input file
    try: param_file
    except:
        stdout.write("parameter file not initiated - check input\n")
        stdout.write(getUsage())
        sys.exit("Incorrect usage!\n")
    try: pqr_file
    except:
        stdout.write("input pqr file not initiated - check input\n")
        stdout.write(getUsage())
        sys.exit("Incorrect usage!\n")
    
    stdout.write("Parsing parameter file %s...\n" % param_file)
    energylist = []
    bradlist = []
    try: f = open(param_file, 'r')
    except IOError:
        stdout.write("cannot open parameter file specified\n")
        sys.exit("IOError!\n")
    index = 0
    for line in f:
        term = line.split()
        if index == 0:
            if term[0].isdigit:
                sdie = float(term[0])
            else:
                stderr.write("main: Parameter format error; First line should be sdie only")
                raise APBSError, "Incorrect format!"
        elif index == 1:
            if len(term) == 1:
                if term[0] == "energy":
                    method = 1
            elif len(term) == 2:
                if term[0] == "radii":
                    column = 0
                else:
                    column = 1
            else:
                stderr.write("main: Parameter format error; Second line should be field flags")
                raise APBSError, "Incorrect format!"
        else:
            if method == 1:
                if term[column].isdigit:
                    energylist.append(float(term[column]))
            else:
                if term[column].isdigit:
                    bradlist.append(float(term[column]))
        index = index + 1
    f.close()
    
    if method:
        numAtoms = len(energylist)
        stdout.write("Parsed energy for %d atoms.\n" % numAtoms)        
    else:
        numAtoms = len(bradlist)
        stdout.write("Parsed Born radii for %d atoms.\n" % numAtoms)
    
    stdout.write("Parsing PQR file %s...\n" % pqr_file)
    position = [[0.0 for i in range(3)] for j in range(numAtoms)]
    x = []
    y = []
    z = []
    chargelist = []
    try: f = open(pqr_file, 'r+')
    except IOError:
        stdout.write("cannot open input pqr file specified\n")
        sys.exit("IOError!\n")
    iatom = 0
    for line in f:
        param = line.split()
        if param[0] != "ATOM":
            continue
        if len(param) == 10:
            x.append(float(param[5]))
            y.append(float(param[6]))
            z.append(float(param[7]))
            chargelist.append(float(param[8])*Python_C)
            iatom=iatom+1
    f.close()
    if len(chargelist) != numAtoms:
        stderr.write("main: Number of atoms in energy list (%d) doest not match PQR list (%d)\n" %
                     (len(chargelist),numAtoms))
        raise APBSError, "Non-matching energy list and PQR file"
    for iatom in range(numAtoms):
        position[iatom][0]=x[iatom]
        position[iatom][1]=y[iatom]
        position[iatom][2]=z[iatom]

    stdout.write("Input files parsed...\n")
    
    # Obtain Born radii from self energies
    stdout.write("Starting energy calculations...\n")
    dij2 = [[0.0 for i in range(numAtoms)] for j in range(numAtoms)]
    fGB = [[0.0 for i in range(numAtoms)] for j in range(numAtoms)]
    if method:
        for i in xrange(numAtoms):
            brad = -pow(chargelist[i],2)*(1-1/sdie)*0.5*Python_Na/(4*pi*Python_e0*energylist[i]*1e3)
            bradlist.append(brad)

    if output_param != "":
        stdout.write("writing parameter file to %s\n" % output_param)
        FILE = open(output_param, "w")
        FILE.write(str(sdie)+"\n")
        if method:
            FILE.write("radii\tenergy\n")
            parameters = zip(bradlist, energylist)
            for i in parameters:
                print >> FILE, "\t".join(map(str,i))
        else:
            FILE.write("radii\n")
            for i in bradlist:
                FILE.write(str(i)+"\n")
        FILE.close()
    
    for i in xrange(numAtoms):
        for j in xrange(i+1):
            
            for coord in xrange(3):
                dij2[i][j] = dij2[i][j] + pow((position[i][coord]-position[j][coord])*1e-10,2)

            d = dij2[i][j]
            bradi = bradlist[i]
            bradj = bradlist[j]
            if j==i:
                fGB[i][j]=bradlist[i]
            else:
                fGB[i][j] = sqrt(d+bradi*bradj*exp(-d/(4.0*bradi*bradj)))
    
    # Calculate energy

    Gpol = 0.0
    for i in xrange(numAtoms):
        for j in xrange(numAtoms):
            if j < i:
                Gpol = Gpol + chargelist[i]*chargelist[j]/fGB[i][j]
            elif j > i:
                Gpol = Gpol + chargelist[i]*chargelist[j]/fGB[j][i]

    for i in xrange(numAtoms):
        Gpol = Gpol + pow(chargelist[i],2)/bradlist[i]
    
    Gpol = -Gpol*(1-1/sdie)*0.5*1e-3*Python_Na/(4*pi*Python_e0)
    
    # Print result
    stdout.write("\nGB Energy: %.10E kJ/mol\n" % Gpol)

    # Record data
    if output_matrix != "":
        FILE = open(filename, 'w')
        term = 0.0
        for i in range(numAtoms):
            for j in range(numAtoms):
                if j<=i:
                    term = chargelist[i]*chargelist[j]/fGB[i][j]
                    term = -term*(1-1/sdie)*0.5*1e-3*Python_Na/(4*pi*Python_e0)
                    FILE.write(str(term)+"\t")
                else:
                    term = chargelist[i]*chargelist[j]/fGB[j][i]
                    term = -term*(1-1/sdie)*0.5*1e-3*Python_Na/(4*pi*Python_e0)
                    FILE.write(str(term)+"\t")
            FILE.write("\n")
        FILE.close()
        stdout.write("Energy matrix output in %s\n" % filename)

    stdout.write("\n")
    stdout.write("Thanks for using APBS!\n\n")

    # Stop the main timer
    main_timer_stop = time.clock()
    stdout.write("Total execution time:  %1.6e sec\n" % (main_timer_stop - main_timer_start))

 
if __name__ == "__main__": main()
