/**
 *  @file    coulomb.c
 *  @author  Nathan Baker
 *  @brief   Small program to calculate Coulombic energies
 *  @version $Id$
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Center for Computational Biology
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2003.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2003.  The Regents of the University of
 * California.  
 * Portions Copyright (c) 1995.  Michael Holst.
 *
 * This file is part of APBS.
 *
 * APBS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * APBS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with APBS; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 * @endverbatim
 */
#include <sys/times.h>
#include <unistd.h>

#include "apbscfg.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vacc.h"

int main(int argc, char **argv) {

    /* OBJECTS */
    Valist *alist;
    Vatom *atom1, *atom2;
    double charge1, charge2, dist, dist2, *pos1, *pos2, energy, myenergy;
    double zmagic, disp[3], force[3], myforce[3];
    int i, j;
    int verbose = 0;
    int doforce = 1;

    /* SYSTEM PARAMETERS */
    char *path;

    char *usage = "\n\n\
This program calculates electrostatic properties using Coulomb's law;\n\
i.e., the Green's function for the Poisson operator in a homogeneous\n\
medium with dielectric constant 1 (vacumm).  It is very important to\n\
realize that all energies, forces, etc. are calculated with this\n\
dielectric constant of 1 and must be scaled to be compared with other\n\
calculations.\n\n\
Usage: coulomb [-v] [-f] <molecule.pqr>\n\n\
   where <molecule.pqr> is the path to the molecule\n\
   structure in PQR format.  This program supports the\n\
   following options:\n\
       -v      give per-atom information\n\
       -f      calculate forces in addition to energies\n\n";

    if ((argc > 4) || (argc < 2)) {
        printf("\n*** Syntax error: got %d arguments, expected 2.\n",argc);
        printf("%s", usage);
        exit(666);
    };
    if (argc > 2) {
        for (i=1; i<(argc-1); i++) {
            if (strcmp("-v", argv[i]) == 0) {
                printf("Providing per-atom information...\n");
                verbose = 1;
            } else if (strcmp("-f", argv[i]) == 0) {
                printf("Calculating forces...\n");
                doforce = 1;
            } else {
                printf("Ignoring option %s\n", argv[i]);
                verbose = 0;
            }
        }
        path = argv[argc-1];
    } else {
        verbose = 0;
        path = argv[1];
    }

    printf("Setting up atom list from %s\n", path);
    alist = Valist_ctor();
    Valist_readPQR(alist, "FILE", "ASC", VNULL, path);
    printf("Read %d atoms\n", Valist_getNumberAtoms(alist));

    /* ENergy scaling factor */
    zmagic  = (1e-3)*(1e10)*Vunit_ec*Vunit_ec*Vunit_Na/(4*VPI*Vunit_eps0);

    energy = 0.0;
    force[0] = 0.0; force[1] = 0.0; force[2] = 0.0;
    
    printf("Using vacuum dielectric, distance in Ang, and charge in e....\n");
    printf("Calculating...\n");
    fflush(stdout);
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom1 = Valist_getAtom(alist, i);
        pos1 = Vatom_getPosition(atom1);
        charge1 = Vatom_getCharge(atom1);
        myenergy = 0.0;
        myforce[0] = 0.0; myforce[1] = 0.0; myforce[2] = 0.0;
        for (j=0; j<Valist_getNumberAtoms(alist); j++) {
            if (j != i) {
                atom2 = Valist_getAtom(alist, j);
                pos2 = Vatom_getPosition(atom2);
                charge2 = Vatom_getCharge(atom2);
                disp[0] = pos1[0] - pos2[0];
                disp[1] = pos1[1] - pos2[1];
                disp[2] = pos1[2] - pos2[2];
                dist2 = (VSQR(disp[0]) + VSQR(disp[1]) + VSQR(disp[2]));
                dist = VSQRT(dist2);
                myenergy += 0.5*(charge1*charge2/dist);
                myforce[0] -= (0.5*disp[0]*(charge1*charge2/(dist*dist2)));
                myforce[1] -= (0.5*disp[1]*(charge1*charge2/(dist*dist2)));
                myforce[2] -= (0.5*disp[2]*(charge1*charge2/(dist*dist2)));
            }
        }
        energy += myenergy;
        force[0] += myforce[0];
        force[1] += myforce[1];
        force[2] += myforce[2];
        myenergy = myenergy*zmagic;
        myforce[0] = myforce[0]*zmagic;
        myforce[1] = myforce[1]*zmagic;
        myforce[2] = myforce[2]*zmagic;
        if (verbose) {
            printf("\tAtom %d:  Energy  = %1.12E kJ/mol\n", i, myenergy);
            if (doforce) {
                printf("\tAtom %d:  x-force = %1.12E kJ/mol/A\n", i, 
                  myforce[0]);
                printf("\tAtom %d:  y-force = %1.12E kJ/mol/A\n", i, 
                  myforce[1]);
                printf("\tAtom %d:  z-force = %1.12E kJ/mol/A\n", i, 
                  myforce[2]);
            }
        }
    }

    energy = energy*zmagic;
    force[0] = force[0]*zmagic;
    force[1] = force[1]*zmagic;
    force[2] = force[2]*zmagic;
 
    printf("\n\n-------------------------------------------------------\n");
    printf("Total energy = %1.12e kJ/mol in vacuum.\n", energy);
    printf("Total x-force = %1.12e kJ/mol/A in vacuum.\n", force[0]);
    printf("Total y-force = %1.12e kJ/mol/A in vacuum.\n", force[1]);
    printf("Total z-force = %1.12e kJ/mol/A in vacuum.\n", force[2]);

    return 0;
}
