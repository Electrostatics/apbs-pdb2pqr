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
 * Copyright (c) 2002-2008, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2008.  Nathan A. Baker
 * Portions Copyright (c) 1999-2002.  The Regents of the University of California.
 * Portions Copyright (c) 1995.  Michael Holst
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  
 * 
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. *
 * @endverbatim
 */

#include "apbscfg.h"
#include "maloc/maloc.h"
#include "apbs/apbs.h"

int main(int argc, char **argv) {

    /* OBJECTS */
    Valist *alist = VNULL;
    Vatom *atom = VNULL;
    Vgreen *green = VNULL;
    Vio *sock = VNULL;

    double *pos, energy, zmagic, disp[3], force[3];
    double *fx, *fy, *fz, *pot, *xp, *yp, *zp, *qp;
    int i, j;
    int doenergy = 0;
    int doforce = 0;

    /* SYSTEM PARAMETERS */
    char *path;

    char *usage = "\n\n\
This program calculates electrostatic properties using Coulomb's law;\n\
i.e., the Green's function for the Poisson operator in a homogeneous\n\
medium with dielectric constant 1 (vacumm).  It is very important to\n\
realize that all energies, forces, etc. are calculated with this\n\
dielectric constant of 1 and must be scaled to be compared with other\n\
calculations.\n\n\
Usage: coulomb [-e] [-f] <molecule.pqr>\n\n\
   where <molecule.pqr> is the path to the molecule\n\
   structure in PQR format.  By default the total energies and forces\n\
   will be printed. This program also supports the following options:\n\
       -e      give per-atom energies\n\
       -f      give per-atom forces\n\n";

    Vio_start();

    if ((argc > 4) || (argc < 2)) {
        Vnm_print(2, "\n*** Syntax error: got %d arguments, expected 2.\n",
           argc);
        Vnm_print(2, "%s", usage);
        exit(666);
    };
    if (argc > 2) {
        for (i=1; i<(argc-1); i++) {
            if (strcmp("-e", argv[i]) == 0) {
                Vnm_print(1, "Providing per-atom energies...\n");
                doenergy = 1;
            } else if (strcmp("-f", argv[i]) == 0) {
                Vnm_print(1, "Providing per-atom forces...\n");
                doforce = 1;
            } else if (strcmp("-ef", argv[i]) == 0 || \
                       strcmp("-fe", argv[i]) == 0){
                Vnm_print(1, "Providing per-atom forces and energies...\n");
                doforce = 1; 
                doenergy = 1;
              
            } else {
                Vnm_print(2, "Ignoring option %s\n", argv[i]);
            }
        }
        path = argv[argc-1];
    } else {
        path = argv[1];
    }

    Vnm_print(1, "Setting up atom list from %s.\n", path);
    alist = Valist_ctor();
    sock = Vio_ctor("FILE", "ASC", VNULL, path, "r");
    if (sock == VNULL) {
        Vnm_print(2, "Problem opening virtual socket %s!\n", 
                  path);
        return 0;
    }
    if (Vio_accept(sock, 0) < 0) {
        Vnm_print(2, "Problem accepting virtual socket %s!\n",
                  path);
        return 0;
    }
    Valist_readPQR(alist,VNULL,sock);
    Vnm_print(1, "Read %d atoms\n", Valist_getNumberAtoms(alist));

    Vnm_print(1, "Setting up Green's function object.\n");
    green = Vgreen_ctor(alist);

    /* Initialize variables */
    Vnm_print(1, "Dielectric constant = 1 (vacuum)\n"); 
    Vnm_print(1, "Distances in Angstroms\n");
    Vnm_print(1, "Charges in electrons\n");
    zmagic  = (1e-3)*Vunit_ec*Vunit_Na;
    energy = 0.0;
    force[0] = 0.0; force[1] = 0.0; force[2] = 0.0;

    Vnm_print(1, "Allocating space for solution...\n");
    fx = Vmem_malloc(VNULL, Valist_getNumberAtoms(alist), sizeof(double));
    fy = Vmem_malloc(VNULL, Valist_getNumberAtoms(alist), sizeof(double));
    fz = Vmem_malloc(VNULL, Valist_getNumberAtoms(alist), sizeof(double));
    pot = Vmem_malloc(VNULL, Valist_getNumberAtoms(alist), sizeof(double));
    xp = Vmem_malloc(VNULL, Valist_getNumberAtoms(alist), sizeof(double));
    yp = Vmem_malloc(VNULL, Valist_getNumberAtoms(alist), sizeof(double));
    zp = Vmem_malloc(VNULL, Valist_getNumberAtoms(alist), sizeof(double));
    qp = Vmem_malloc(VNULL, Valist_getNumberAtoms(alist), sizeof(double));
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        fx[i] = 0.0; 
        fy[i] = 0.0; 
        fz[i] = 0.0;
        pot[i] = 0.0;
        atom = Valist_getAtom(alist, i);
        pos = Vatom_getPosition(atom);
        xp[i] = pos[0]; 
        yp[i] = pos[1]; 
        zp[i] = pos[2];
        qp[i] = Vatom_getCharge(atom);
    }

    Vnm_print(1, "Calculating...\n");
    Vgreen_coulombD(green, Valist_getNumberAtoms(alist), xp, yp, zp, 
            pot, fx, fy, fz);

    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        fx[i] *= (-0.5*qp[i]*zmagic);
        fy[i] *= (-0.5*qp[i]*zmagic);
        fz[i] *= (-0.5*qp[i]*zmagic);
        pot[i] *= (0.5*qp[i]*zmagic);
        energy += pot[i];
        force[0] += fx[i];
        force[1] += fy[i];
        force[2] += fz[i];
        if (doenergy) {
            Vnm_print(1, "\tAtom %d:  Energy  = %1.12E kJ/mol\n", i+1, pot[i]);
        }
        if (doforce) {
          Vnm_print(1, "\tAtom %d:  x-force = %1.12E kJ/mol/A\n", i+1, 
                    fx[i]);
          Vnm_print(1, "\tAtom %d:  y-force = %1.12E kJ/mol/A\n", i+1, 
                    fy[i]);
          Vnm_print(1, "\tAtom %d:  z-force = %1.12E kJ/mol/A\n", i+1, 
                    fz[i]);
        }    
    }

    Vnm_print(1, "\n\n-------------------------------------------------------\n");
    Vnm_print(1, "Total energy = %1.12e kJ/mol in vacuum.\n", energy);
    Vnm_print(1, "Total x-force = %1.12e kJ/mol/A in vacuum.\n", force[0]);
    Vnm_print(1, "Total y-force = %1.12e kJ/mol/A in vacuum.\n", force[1]);
    Vnm_print(1, "Total z-force = %1.12e kJ/mol/A in vacuum.\n", force[2]);

    Vmem_free(VNULL, Valist_getNumberAtoms(alist), sizeof(double), 
      (void **)&fx);
    Vmem_free(VNULL, Valist_getNumberAtoms(alist), sizeof(double), 
      (void **)&fy);
    Vmem_free(VNULL, Valist_getNumberAtoms(alist), sizeof(double), 
      (void **)&fz);
    Vmem_free(VNULL, Valist_getNumberAtoms(alist), sizeof(double), 
      (void **)&pot);
    Vmem_free(VNULL, Valist_getNumberAtoms(alist), sizeof(double), 
      (void **)&xp);
    Vmem_free(VNULL, Valist_getNumberAtoms(alist), sizeof(double), 
      (void **)&yp);
    Vmem_free(VNULL, Valist_getNumberAtoms(alist), sizeof(double), 
      (void **)&zp);
    Vmem_free(VNULL, Valist_getNumberAtoms(alist), sizeof(double), 
      (void **)&qp);
    Vgreen_dtor(&green);
    Valist_dtor(&alist);

    return 0;
}
