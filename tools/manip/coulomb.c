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
 * Nathan A. Baker (nbaker@wasabi.ucsd.edu)
 * Dept. of Chemistry and Biochemistry
 * University of California, San Diego 
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 1999-2002.  Nathan A. Baker.  All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above copyright notice, this paragraph and the
 * following two paragraphs appear in all copies, modifications, and
 * distributions.
 *
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE
 * AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE AUTHORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
 * HEREUNDER IS PROVIDED "AS IS".  THE AUTHORS HAVE NO OBLIGATION TO PROVIDE
 * MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
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
    double charge1, charge2, dist, *pos1, *pos2, energy;
    int i, j;

    /* SYSTEM PARAMETERS */
    char *path;

    char *usage = "\n\nUsage: coulomb <molecule.pqr>\n\n"
                  "\twhere <molecule.pqr> is the path to the molecule structure in PQR format\n\n";

    if (argc != 2) {
        printf("\n*** Syntax error: got %d arguments, expected 2.\n",argc);
        printf("%s", usage);
        exit(666);
    };
    path = argv[1];

    printf("Setting up atom list from %s\n", path);
    alist = Valist_ctor();
    Valist_readPQR(alist, "FILE", "ASC", VNULL, path);
    printf("Read %d atoms\n", Valist_getNumberAtoms(alist));

    energy = 0.0;

    printf("Calculating energy of system in vacuum (distances in Ang)...\n");
    fflush(stdout);
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom1 = Valist_getAtom(alist, i);
        pos1 = Vatom_getPosition(atom1);
        charge1 = Vatom_getCharge(atom1);
        for (j=i+1; j<Valist_getNumberAtoms(alist); j++) {
            atom2 = Valist_getAtom(alist, j);
            pos2 = Vatom_getPosition(atom2);
            charge2 = Vatom_getCharge(atom2);
            dist = VSQRT(VSQR(pos1[0]-pos2[0]) + VSQR(pos1[1]-pos2[1])
                               + VSQR(pos1[2]-pos2[2]));
            energy += (charge1*charge2/dist);
        }
    }

    energy = energy*(1e-3)*(1e10)*Vunit_ec*Vunit_ec*Vunit_Na/(4*VPI*Vunit_eps0);
 
    printf("Energy = %1.12e kJ/mol in vacuum.\n", energy);

    return 0;
}
