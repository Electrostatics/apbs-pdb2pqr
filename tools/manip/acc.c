/**
 *  @file    acc.c
 *  @author  Nathan Baker
 *  @brief   Small program to calculate volumes, areas, etc. of molecules
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
 * Copyright (c) 2002-2004.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2002.  The Regents of the University of
 * California.  
 * Portions Copyright (c) 1995.  Michael Holst.
 *
 * This file is part of APBS.
 *
 * APBS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public Liorise as published by
 * the Free Software Foundation; either version 2 of the Liorise, or
 * (at your option) any later version.
 *
 * APBS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public Liorise for more details.
 *
 * You should have received a copy of the GNU General Public Liorise
 * along with APBS; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 * @endverbatim
 */

#include "apbscfg.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vclist.h"
#include "apbs/vacc.h"

/** 
 * @brief  Print usage information
 */
int usage(
        int rc /** Return code, passed through */
        ) {

    char *ustr = \
      "\nUsage:  acc [options] <molecule.pqr>\n\n"
      "\tThe [options] arguments can include any of the following:\n"
      "\t\t--probe=<value>  Specify the probe radius (in Angstroms).\n"
      "\t\t  Default = 1.4 A.\n"
      "\t\t--vol-density=<value>  Specify the density of grid points to\n"
      "\t\t  use for volume quadratures.  Default = 2.0 A^{-3}.\n"
      "\t\t--surf-density=<value>  Specify the density of grid points to\n"
      "\t\t  use for area quadratures.  Default = 1.0 A^{-2}\n"
      "\t\t--verbose  Increase the verbosity of the output to include\n"
      "\t\t  per-atom information, where applicable\n"
      "\t\t--area-only  Only calculate the surface areas\n"
      "\t\t--help,-h  Print this message\n"
      "\tand <molecule.pqr> is the path to the molecule structure in PQR\n"
      "\tformat\n\n";

    Vnm_print(2, "%s", ustr);

    return rc;
}

/**
 * @brief Main code
 */
int main(int argc, char **argv) {

    /* OBJECTS */
    Valist *alist;
    Vclist *clist;
    Vacc  *acc;
    Vatom *atom;
    Vio *sock = VNULL;
    char *substr;
    char *path;

    /* VCLIST PARAMETERS */
    int nhash[3] = {60, 60, 60};

    /* QUADRATURE VARIABLES */
    double vdwVol = 0.0;
    double ivdwVol = 0.0;
    double molVol = 0.0;
    double sasa = 0.0;
    double atom_sasa = 0.0;

    /* Quadrature steps */
    int i, npts[3];
    double spacs[3], vec[3];
    double w, wx, wy, wz, len, fn, x, y, z, vol;
    double *lower_corner, *upper_corner;

    /* Default parameters */
    double probe_radius = 1.4;
    double vol_density = 2.0;
    double surf_density = 1.0;
    double fVerbose = 0;
    int fAreaOnly = 0;

    Vio_start();

    /* Check usage */
    if (argc == 1) {
        Vnm_print(2, "\nError:  didn't get any arguments!\n");
        return usage(EXIT_FAILURE);
    }
    /* Check for help */
    for (i=1; i<argc; i++) {
        if (strstr(argv[i], "--help") != NULL) return usage(EXIT_SUCCESS);
        if (strstr(argv[i], "-h") != NULL) return usage(EXIT_SUCCESS);
    }
    /* Molecule path */
    path = argv[argc-1];
    /* Remaining options */
    for (i=1; i<(argc-1); i++) {
        /* Look for probe radius specification */
        if (strstr(argv[i], "--probe=") != NULL) {
            substr = strchr(argv[i], '=');
            substr = substr + 1;
            if (sscanf(substr, "%lf", &probe_radius) == 0) {
                Vnm_print(2, 
                        "\nError:  unable to parse (%s) as float!\n", 
                        substr); 
                return usage(EXIT_FAILURE);
            } 
        } else if (strstr(argv[i], "--vol-density=") != NULL) {
            substr = strchr(argv[i], '=');
            substr = substr + 1;
            if (sscanf(substr, "%lf", &vol_density) == 0) {
                Vnm_print(2, 
                        "\nError:  unable to parse (%s) as float!\n", 
                        substr); 
                return usage(EXIT_FAILURE);
            } 
        } else if (strstr(argv[i], "--surf-density=") != NULL) {
            substr = strchr(argv[i], '=');
            substr = substr + 1;
            if (sscanf(substr, "%lf", &surf_density) == 0) {
                Vnm_print(2, 
                        "\nError:  unable to parse (%s) as float!\n", 
                        substr); 
                return usage(EXIT_FAILURE);
            } 
        } else if (strstr(argv[i], "--verbose") != NULL) {
            fVerbose = 1;
        } else if (strstr(argv[i], "--area-only") != NULL) {
            fAreaOnly = 1;
        } else {
            Vnm_print(2, "\nError:  unknown option (%s)\n", argv[i]);
            return usage(EXIT_FAILURE);
        }
    }

    /* Parameters */
    Vnm_print(1, "Molecule path:  %s\n", path);
    Vnm_print(1, "Probe radius:  %g A\n", probe_radius);
    Vnm_print(1, "Volume point density:  %g A^{-3}\n", vol_density);
    Vnm_print(1, "Surface point density:  %g A^{-2}\n", surf_density);
    if (fVerbose) Vnm_print(1, "Verbose output.\n");
    if (fAreaOnly) Vnm_print(1, "Area-only output.\n");

    /* Read atom list */
    Vnm_print(1, "\nReading PQR file...\n");
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
    Valist_readPQR(alist,sock);

    /* Set up Vacc and Vclist */
    Vnm_print(1, "Setting up hash table and accessibility object...\n");
    clist = Vclist_ctor(alist, probe_radius, nhash, CLIST_AUTO_DOMAIN, 
            VNULL, VNULL);
    acc = Vacc_ctor(alist, clist, surf_density);

    if (!fAreaOnly) {

        /* Set up quadrature */
        lower_corner = clist->lower_corner;
        upper_corner = clist->upper_corner;
        vol = 1.0;
        for (i=0; i<3; i++) {
            len = upper_corner[i] - lower_corner[i];
            vol *= len;
            fn = len*vol_density + 1;
            npts[i] = (int)ceil(fn);
            spacs[i] = len/((double)(npts[i])-1.0);
        }

        Vnm_print(1, "Quadrature mesh spacing = (%g, %g, %g)\n",
                spacs[0], spacs[1], spacs[2]);
        Vnm_print(1, "Quadrature mesh points = (%d, %d, %d)\n",
                npts[0], npts[1], npts[2]);

        Vnm_print(1, "\nPerforming volume quadrature...\n");

        for (x=lower_corner[0]; x<=upper_corner[0]; x=x+spacs[0]) {
            if ( VABS(x - lower_corner[0]) < VSMALL) {
                wx = 0.5;
            } else if ( VABS(x - upper_corner[0]) < VSMALL) {
                wx = 0.5;
            } else {
                wx = 1.0;
            }
            vec[0] = x;
            for (y=lower_corner[1]; y<=upper_corner[1]; y=y+spacs[1]) {
                if ( VABS(y - lower_corner[1]) < VSMALL) {
                    wy = 0.5;
                } else if ( VABS(y - upper_corner[1]) < VSMALL) {
                    wy = 0.5;
                } else {
                    wy = 1.0;
                }
                vec[1] = y;
                for (z=lower_corner[2]; z<=upper_corner[2]; z=z+spacs[2]) {
                    if ( VABS(z - lower_corner[2]) < VSMALL) {
                        wz = 0.5;
                    } else if ( VABS(z - upper_corner[2]) < VSMALL) {
                        wz = 0.5;
                    } else {
                        wz = 1.0;
                    }
                    vec[2] = z;

                    w = wx*wy*wz;
                    
                    /* printf("%g, %g, %g (%g)\n", x, y, z, w); */

                    vdwVol += (w*(1.0-Vacc_vdwAcc(acc, vec)));
                    ivdwVol += (w*(1.0-Vacc_ivdwAcc(acc, vec, probe_radius)));
                    molVol += (w*(1.0-Vacc_molAcc(acc, vec, probe_radius)));

                } /* z loop */
            } /* y loop */
        } /* x loop */

        w  = spacs[0]*spacs[1]*spacs[2];
        vdwVol *= w;
        ivdwVol *= w;
        molVol *= w;

        Vnm_print(1, "van der Waals volume = %g A^3\n", vdwVol);
        Vnm_print(1, "Inflated van der Waals volume = %g A^3\n", ivdwVol);
        Vnm_print(1, "Molecular volume = %g A^3\n", molVol);

    } /* if !fAreaOnly */

    Vnm_print(1, "\nCalculating solvent accessible surface areas...\n");
    if (fVerbose) {

        sasa = 0.0;
        Vnm_print(1, "Atom\tArea (A^2)\n");
        Vnm_print(1, "----\t------------------\n");
        for (i=0; i<Valist_getNumberAtoms(alist); i++) {
            atom = Valist_getAtom(alist, i);
            atom_sasa = Vacc_atomSASA(acc, probe_radius, atom);
            sasa += atom_sasa;
            Vnm_print(1, "%d\t%1.12E\n", atom->id, atom_sasa);
        }
        Vnm_print(1, "----\t------------------\n");
        Vnm_print(1, "TOTL\t%1.12E\n", sasa);
    } else {
        Vnm_print(1, "Total SASA:  %1.12E\n", 
                Vacc_totalSASA(acc, probe_radius));
    } /* if fVerbose */

    return EXIT_SUCCESS;
}
