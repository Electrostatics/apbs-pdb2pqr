/**
 *  @file    routines.c
 *  @author  Nathan Baker
 *  @brief   Supporting routines for APBS front end
 *  @version $Id$
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2003.  Washington University in St. Louis.
 * All Rights Reserved.
 *
 * Portions Copyright (c) 1999-2003.  The Regents of the University of
 * California.  
 * Portions Copyright (c) 1995.  Michael Holst.
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

#include "apbscfg.h"
#include "apbs/apbs.h"  
#include "apbs/nosh.h"  
#include "apbs/vgrid.h"  
#include "apbs/mgparm.h"  
#include "apbs/pbeparm.h"  
#include "apbs/femparm.h"  
#include "maloc/maloc.h"  

#include "routines.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  startVio
//
// Purpose:  Wrapper for Vio_start
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void startVio() { Vio_start(); }


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  loadMolecules
//
// Purpose:  Load molecules from files
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int loadMolecules(NOsh *nosh, Valist *alist[NOSH_MAXMOL]) {
    
    int i, j;
    double q; 
    Vatom *atom;

    Vnm_tprint( 1, "Got PQR paths for %d molecules\n", nosh->nmol);
    if (nosh->nmol <= 0) {
       Vnm_tprint(2, "You didn't specify any molecules (correctly)!\n");
       Vnm_tprint(2, "Bailing out!\n");
       return 0;
    }

    for (i=0; i<nosh->nmol; i++) {
        Vnm_tprint( 1, "Reading atom data from %s:\n",
          nosh->molpath[i]);
        alist[i] = Valist_ctor();
        if (Valist_readPQR(alist[i], "FILE", "ASC", VNULL,
          nosh->molpath[i]) != 1) {
            Vnm_tprint( 2, "Fatal error while reading from %s\n",
              nosh->molpath[i]);
            return 0;
        } else {
            Vnm_tprint( 1, "  %d atoms\n",
              Valist_getNumberAtoms(alist[i]));
            Vnm_tprint( 1, "  Centered at (%4.3e, %4.3e, %4.3e)\n",
              alist[i]->center[0], alist[i]->center[1], alist[i]->center[2]);
            Vnm_tprint( 1, "  Net charge %3.2e e\n",
              alist[i]->charge);        
            /* Check for uncharged molecule */
            q = 0;
            for (j=0; j<Valist_getNumberAtoms(alist[i]); j++) {
                atom = Valist_getAtom(alist[i], j);
                q += VSQR(Vatom_getCharge(atom));
            }
            if (q < (1e-6)) {
                Vnm_print(2, "Molecule #%d is uncharged!\n");
                Vnm_print(2, "Sum square charge = %g\n", q);
                Vnm_print(2, "Bailing out!\n");
            }
        }
    }

    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  killMolecules
//
// Purpose:  Kill the molecule structures
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void killMolecules(NOsh *nosh, Valist *alist[NOSH_MAXMOL]) {
    
    int i;

#ifndef VAPBSQUIET
    Vnm_tprint( 1, "Destroying %d molecules\n", nosh->nmol);
#endif

    for (i=0; i<nosh->nmol; i++) Valist_dtor(&(alist[i]));

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  loadDielMaps
// 
// Purpose:  Load the dielectric maps from files
// 
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int loadDielMaps(NOsh *nosh, 
  Vgrid *dielXMap[NOSH_MAXMOL], Vgrid *dielYMap[NOSH_MAXMOL],
  Vgrid *dielZMap[NOSH_MAXMOL]) {

    int i, ii, nx, ny, nz;
    double sum, hx, hy, hzed, xmin, ymin, zmin;

    if (nosh->ndiel > 0) 
      Vnm_tprint( 1, "Got paths for %d dielectric map sets\n", 
        nosh->ndiel);
    else return 1;

    for (i=0; i<nosh->ndiel; i++) {
        Vnm_tprint( 1, "Reading x-shifted dielectric map data from \
%s:\n", nosh->dielXpath[i]);
        dielXMap[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        if (nosh->dielfmt[i] == 0) {
            if (Vgrid_readDX(dielXMap[i], "FILE", "ASC", VNULL, 
              nosh->dielXpath[i]) != 1) {
                Vnm_tprint( 2, "Fatal error while reading from %s\n",
                  nosh->dielXpath[i]);
                return 0;
            }
            nx = dielXMap[i]->nx;
            ny = dielXMap[i]->ny;
            nz = dielXMap[i]->nz;
            hx = dielXMap[i]->hx;
            hy = dielXMap[i]->hy;
            hzed = dielXMap[i]->hzed;
            xmin = dielXMap[i]->xmin;
            ymin = dielXMap[i]->ymin;
            zmin = dielXMap[i]->zmin;
            Vnm_tprint(1, "  %d x %d x %d grid\n", nx, ny, nz);
            Vnm_tprint(1, "  (%g, %g, %g) A spacings\n", hx, hy, hzed);
            Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n", 
              xmin, ymin, zmin);
            sum = 0;
            for (ii=0; ii<(nx*ny*nz); ii++)
              sum += (dielXMap[i]->data[ii]);
            sum = sum*hx*hy*hzed;
            Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", 
              sum);
        } else {
            Vnm_tprint( 2, "INVALID FORMAT!\n");
            return 0;
        }
        Vnm_tprint( 1, "Reading y-shifted dielectric map data from \
%s:\n", nosh->dielYpath[i]);
        dielYMap[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        if (nosh->dielfmt[i] == 0) {
            if (Vgrid_readDX(dielYMap[i], "FILE", "ASC", VNULL, nosh->dielYpath[i])
              != 1) {
                Vnm_tprint( 2, "Fatal error while reading from %s\n",
                  nosh->dielYpath[i]);
                return 0;
            }
            nx = dielYMap[i]->nx;
            ny = dielYMap[i]->ny;
            nz = dielYMap[i]->nz;
            hx = dielYMap[i]->hx;
            hy = dielYMap[i]->hy;
            hzed = dielYMap[i]->hzed;
            xmin = dielYMap[i]->xmin;
            ymin = dielYMap[i]->ymin;
            zmin = dielYMap[i]->zmin;
            Vnm_tprint(1, "  %d x %d x %d grid\n", nx, ny, nz);
            Vnm_tprint(1, "  (%g, %g, %g) A spacings\n", hx, hy, hzed);
            Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
              xmin, ymin, zmin);
            sum = 0;
            for (ii=0; ii<(nx*ny*nz); ii++)
              sum += (dielYMap[i]->data[ii]);
            sum = sum*hx*hy*hzed;
            Vnm_tprint(1, "  Volume integral = %3.2e A^3\n",
              sum);
        } else {
            Vnm_tprint( 2, "INVALID FORMAT!\n");
            return 0;
        }
        Vnm_tprint( 1, "Reading z-shifted dielectric map data from \
%s:\n", nosh->dielZpath[i]);
        dielZMap[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        if (nosh->dielfmt[i] == 0) {
            if (Vgrid_readDX(dielZMap[i], "FILE", "ASC", VNULL, nosh->dielZpath[i])
              != 1) {
                Vnm_tprint( 2, "Fatal error while reading from %s\n",
                  nosh->dielZpath[i]);
                return 0;
            }
            nx = dielZMap[i]->nx;
            ny = dielZMap[i]->ny;
            nz = dielZMap[i]->nz;
            hx = dielZMap[i]->hx;
            hy = dielZMap[i]->hy;
            hzed = dielZMap[i]->hzed;
            xmin = dielZMap[i]->xmin;
            ymin = dielZMap[i]->ymin;
            zmin = dielZMap[i]->zmin;
            Vnm_tprint(1, "  %d x %d x %d grid\n",
              nx, ny, nz);
            Vnm_tprint(1, "  (%g, %g, %g) A spacings\n",
              hx, hy, hzed);
            Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
              xmin, ymin, zmin);
            sum = 0;
            for (ii=0; ii<(nx*ny*nz); ii++)
              sum += (dielZMap[i]->data[ii]);
            sum = sum*hx*hy*hzed;
            Vnm_tprint(1, "  Volume integral = %3.2e A^3\n",
              sum);
        } else {
            Vnm_tprint( 2, "INVALID FORMAT!\n");
            return 0;
        }
    }

    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  killDielMaps
// 
// Purpose:  Kill the dielectric map structures
// 
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void killDielMaps(NOsh *nosh, 
  Vgrid *dielXMap[NOSH_MAXMOL], Vgrid *dielYMap[NOSH_MAXMOL],
  Vgrid *dielZMap[NOSH_MAXMOL]) {

    int i;

    if (nosh->ndiel > 0) {
#ifndef VAPBSQUIET
	Vnm_tprint( 1, "Destroying %d dielectric map sets\n", 
		    nosh->ndiel);
#endif
	for (i=0; i<nosh->ndiel; i++) {
	    Vgrid_dtor(&(dielXMap[i]));
	    Vgrid_dtor(&(dielYMap[i]));
	    Vgrid_dtor(&(dielZMap[i]));
	}
    }
    else return;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  loadKappaMaps
//
// Purpose:  Load kappa maps from files
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int loadKappaMaps(NOsh *nosh, Vgrid *map[NOSH_MAXMOL]) {

    int i, ii;
    double sum;

    if (nosh->nkappa > 0) 
      Vnm_tprint( 1, "Got paths for %d kappa maps\n", nosh->nkappa);
    else return 1;

    for (i=0; i<nosh->nkappa; i++) {
        Vnm_tprint( 1, "Reading kappa map data from %s:\n",
          nosh->kappapath[i]);
        map[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        if (nosh->kappafmt[i] == 0) {
            if (Vgrid_readDX(map[i], "FILE", "ASC", VNULL, nosh->kappapath[i]) 
              != 1) {
                Vnm_tprint( 2, "Fatal error while reading from %s\n",
                  nosh->kappapath[i]);
                return 0;
            }
            Vnm_tprint(1, "  %d x %d x %d grid\n", 
              map[i]->nx, map[i]->ny, map[i]->nz);
            Vnm_tprint(1, "  (%g, %g, %g) A spacings\n", 
              map[i]->hx, map[i]->hy, map[i]->hzed);
            Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n", 
              map[i]->xmin, map[i]->ymin, map[i]->zmin);
            sum = 0;
            for (ii=0; ii<(map[i]->nx*map[i]->ny*map[i]->nz); ii++)
              sum += (map[i]->data[ii]);
            sum = sum*map[i]->hx*map[i]->hy*map[i]->hzed;
            Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);

        } else {
            Vnm_tprint( 2, "INVALID FORMAT!\n");
            return 0;
        }
    }

    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  killKappaMaps
//
// Purpose:  Kill kappa map structures
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void killKappaMaps(NOsh *nosh, Vgrid *map[NOSH_MAXMOL]) {

    int i;

    if (nosh->nkappa > 0) {
#ifndef VAPBSQUIET
      Vnm_tprint( 1, "Destroying %d kappa maps\n", nosh->nkappa);
#endif
      for (i=0; i<nosh->nkappa; i++) Vgrid_dtor(&(map[i]));
    }
    else return;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  loadChargeMaps
// 
// Purpose:  Load charge maps from files
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int loadChargeMaps(NOsh *nosh, Vgrid *map[NOSH_MAXMOL]) {

    int i, ii;
    double sum;

    if (nosh->ncharge > 0)
      Vnm_tprint( 1, "Got paths for %d charge maps\n", nosh->ncharge);
    else return 1;

    for (i=0; i<nosh->ncharge; i++) {
        Vnm_tprint( 1, "Reading charge map data from %s:\n",
          nosh->chargepath[i]);
        map[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        if (nosh->chargefmt[i] == 0) {
            if (Vgrid_readDX(map[i], "FILE", "ASC", VNULL, nosh->chargepath[i])
              != 1) {
                Vnm_tprint( 2, "Fatal error while reading from %s\n",
                  nosh->chargepath[i]);
                return 0;
            }
            Vnm_tprint(1, "  %d x %d x %d grid\n", 
              map[i]->nx, map[i]->ny, map[i]->nz);
            Vnm_tprint(1, "  (%g, %g, %g) A spacings\n", 
              map[i]->hx, map[i]->hy, map[i]->hzed);
            Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n", 
              map[i]->xmin, map[i]->ymin, map[i]->zmin);
            sum = 0;
            for (ii=0; ii<(map[i]->nx*map[i]->ny*map[i]->nz); ii++) 
              sum += (map[i]->data[ii]);
            sum = sum*map[i]->hx*map[i]->hy*map[i]->hzed;
            Vnm_tprint(1, "  Charge map integral = %3.2e e\n", sum);
        } else {
            Vnm_tprint( 2, "INVALID FORMAT!\n");
            return 0;
        }
    }

    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  killChargeMaps
// 
// Purpose:  Kill charge map structures
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void killChargeMaps(NOsh *nosh, Vgrid *map[NOSH_MAXMOL]) {

    int i;

    if (nosh->ncharge > 0) {
#ifndef VAPBSQUIET
      Vnm_tprint( 1, "Destroying %d charge maps\n", nosh->ncharge);
#endif

      for (i=0; i<nosh->ncharge; i++) Vgrid_dtor(&(map[i]));
    }

    else return;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  printPBEPARM
//
// Purpose:  Print useful stuff from the PBE parameter file
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void printPBEPARM(PBEparm *pbeparm) {
    
    int i;
    double ionstr = 0.0;

    for (i=0; i<pbeparm->nion; i++)
      ionstr += 0.5*(VSQR(pbeparm->ionq[i])*pbeparm->ionc[i]);

    Vnm_tprint( 1, "  Molecule ID: %d\n", pbeparm->molid);
    if (pbeparm->nonlin) Vnm_tprint( 1, "  Nonlinear PBE\n");
    else Vnm_tprint( 1, "  Linearized PBE\n");
    if (pbeparm->bcfl == 0) {
        Vnm_tprint( 1, "  Zero boundary conditions\n");
    } else if (pbeparm->bcfl == 1) {
        Vnm_tprint( 1, "  Single Debye-Huckel sphere boundary \
conditions\n");
    } else if (pbeparm->bcfl == 2) {
        Vnm_tprint( 1, "  Multiple Debye-Huckel sphere boundary \
conditions\n");
    } else if (pbeparm->bcfl == 4) {
        Vnm_tprint( 1, "  Boundary conditions from focusing\n");
    }
    Vnm_tprint( 1, "  %d ion species (%4.3f M ionic strength):\n",
      pbeparm->nion, ionstr);
    for (i=0; i<pbeparm->nion; i++) {
        Vnm_tprint( 1, "    %4.3f A-radius, %4.3f e-charge, \
%4.3f M concentration\n", 
          pbeparm->ionr[i], pbeparm->ionq[i], pbeparm->ionc[i]);            
    }
    Vnm_tprint( 1, "  Solute dielectric: %4.3f\n", pbeparm->pdie);
    Vnm_tprint( 1, "  Solvent dielectric: %4.3f\n", pbeparm->sdie);
    switch (pbeparm->srfm) {
        case 0:
            Vnm_tprint( 1, "  Using \"molecular\" surface \
definition; no smoothing\n");
            Vnm_tprint( 1, "  Solvent probe radius: %4.3f A\n",
              pbeparm->srad);
            break;
        case 1:
            Vnm_tprint( 1, "  Using \"molecular\" surface definition;\
 harmonic average smoothing\n");
            Vnm_tprint( 1, "  Solvent probe radius: %4.3f A\n",
              pbeparm->srad);
            break;
        case 2:
            Vnm_tprint( 1, "  Using spline-based surface definition;\
 window = %4.3f\n", pbeparm->swin);
            break;
        default:
            break;
    }
    switch (pbeparm->chgm) {
        case 0:
            Vnm_tprint(1, "  Using linear spline charge discretization.\n");
            break;
        case 1:
            Vnm_tprint(1, "  Using cubic spline charge discretization.\n");
            break;
        default:
            break;
    }
    Vnm_tprint( 1, "  Temperature:  %4.3f K\n", pbeparm->temp);
    Vnm_tprint( 1, "  Surface tension:  %4.3f kJ/mol/A^2\n",
      pbeparm->gamma);
    if (pbeparm->calcenergy == 1) Vnm_tprint( 1, "  Electrostatic \
energies will be calculated\n");
    if (pbeparm->calcforce == 1) Vnm_tprint( 1, "  Net solvent \
forces will be calculated \n");
    if (pbeparm->calcforce == 2) Vnm_tprint( 1, "  All-atom \
solvent forces will be calculated\n");
    for (i=0; i<pbeparm->numwrite; i++) {
        switch (pbeparm->writetype[i]) {
            case VDT_CHARGE:
                Vnm_tprint(1, "  Charge distribution to be written to ");
                break;
            case VDT_POT:
                Vnm_tprint(1, "  Potential to be written to ");
                break;
            case VDT_SMOL:
                Vnm_tprint(1, "  Molecular solvent accessibility \
to be written to ");
                break;
            case VDT_SSPL:
                Vnm_tprint(1, "  Spline-based solvent accessibility \
to be written to ");
                break;
            case VDT_VDW:
                Vnm_tprint(1, "  van der Waals solvent accessibility \
to be written to ");
                break;
            case VDT_IVDW:
                Vnm_tprint(1, "  Ion accessibility to be written to ");
                break;
            case VDT_LAP:
                Vnm_tprint(1, "  Potential Laplacian to be written to ");
                break;
            case VDT_EDENS:
                Vnm_tprint(1, "  Energy density to be written to ");
                break;
            case VDT_NDENS:
                Vnm_tprint(1, "  Ion number density to be written to ");
                break;
            case VDT_QDENS:
                Vnm_tprint(1, "  Ion charge density to be written to ");
                break;
            case VDT_DIELX:
                Vnm_tprint(1, "  X-shifted dielectric map to be written \
to ");
                break;
            case VDT_DIELY:
                Vnm_tprint(1, "  Y-shifted dielectric map to be written \
to ");
                break;
            case VDT_DIELZ:
                Vnm_tprint(1, "  Z-shifted dielectric map to be written \
to ");
                break;
            case VDT_KAPPA:
                Vnm_tprint(1, "  Kappa map to be written to ");
                break;
            default: 
                Vnm_tprint(2, "  Invalid data type for writing!\n");
                break;
        }
        switch (pbeparm->writefmt[i]) {
            case VDF_DX:
                Vnm_tprint(1, "%s.%s\n", pbeparm->writestem[i], "dx");
                break;
            case VDF_UHBD:
                Vnm_tprint(1, "%s.%s\n", pbeparm->writestem[i], "grd");
                break;
            case VDF_AVS:
                Vnm_tprint(1, "%s.%s\n", pbeparm->writestem[i], "ucd");
                break;
            default: 
                Vnm_tprint(2, "  Invalid format for writing!\n");
                break;
        }
 
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  printMGPARM
//
// Purpose:  Print useful stuff from the MG parameter file
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void printMGPARM(MGparm *mgparm, double realCenter[3]) {

    if (mgparm->type == 2) {
        Vnm_tprint( 1, "  Partition overlap fraction = %g\n", 
          mgparm->ofrac);
        Vnm_tprint( 1, "  Processor array = %d x %d x %d\n", 
          mgparm->pdime[0], mgparm->pdime[1], mgparm->pdime[2]);
    }
    Vnm_tprint( 1, "  Grid dimensions: %d x %d x %d\n",
      mgparm->dime[0], mgparm->dime[1], mgparm->dime[2]);
    Vnm_tprint( 1, "  Grid spacings: %4.3f x %4.3f x %4.3f\n",
      mgparm->grid[0], mgparm->grid[1], mgparm->grid[2]);
    Vnm_tprint( 1, "  Grid lengths: %4.3f x %4.3f x %4.3f\n",
      mgparm->glen[0], mgparm->glen[1], mgparm->glen[2]);
    Vnm_tprint( 1, "  Grid center: (%4.3f, %4.3f, %4.3f)\n",
      realCenter[0], realCenter[1], realCenter[2]);
    Vnm_tprint( 1, "  Multigrid levels: %d\n", mgparm->nlev);

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  initMG
//
// Purpose:  Setup a MG calculation
//
// Args:     realCenter    The actual center of the fine mesh (this could be
//                         somewhat different than a molecule center in the
//                         case of parallel focusing)
//           nosh          Holds input file
//           pbeparm       PBE parameters for this calc
//           mgparm        Multigrid parameters for this calc
//           pbe           PBE object (accessibility, etc. inside)
//           pmgp          Array of PMG parameter objects
//           pmg           Array of PMG objects
//           i             Index of this calculation in pmgp/pmg arrays
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int initMG(int i, NOsh *nosh, MGparm *mgparm, 
  PBEparm *pbeparm, double realCenter[3], Vpbe *pbe[NOSH_MAXCALC], 
  Valist *alist[NOSH_MAXMOL], Vgrid *dielXMap[NOSH_MAXMOL], 
  Vgrid *dielYMap[NOSH_MAXMOL], Vgrid *dielZMap[NOSH_MAXMOL],
  Vgrid *kappaMap[NOSH_MAXMOL], Vgrid *chargeMap[NOSH_MAXMOL],
  Vpmgp *pmgp[NOSH_MAXCALC], Vpmg *pmg[NOSH_MAXCALC]) {
    
    int j, bytesTotal, highWater, imol;
    double sparm, iparm;
    Vgrid *theDielXMap, *theDielYMap, *theDielZMap, *theKappaMap, *theChargeMap;

    Vnm_tstart(27, "Setup timer");

    /* Fix mesh center for "GCENT MOL #" types of declarations. */
    if (mgparm->cmeth == 1) {
        for (j=0; j<3; j++) {
            imol = mgparm->centmol-1;
            if (imol < nosh->nmol) {
                mgparm->center[j] = (alist[imol])->center[j];
            } else{ 
                Vnm_tprint(2, "ERROR!  Bogus molecule number (%d) for \
fgcent/cgcent!\n",  (imol+1));
                return 0;
            }
        }
    }

    /* If we're a parallel calculation, update the grid center based on
     * the appropriate shifts */
    if (mgparm->type == 2) {
        for (j=0; j<3; j++) realCenter[j] = mgparm->center[j]
          + mgparm->partOlapCenterShift[j];
    } else {
        for (j=0; j<3; j++) realCenter[j] = mgparm->center[j];
    }

    /* Set up PBE object */
    if (pbeparm->srfm == 2) sparm = pbeparm->swin;
    else sparm = pbeparm->srad;
    if (pbeparm->nion > 0) iparm = pbeparm->ionr[0];
    else iparm = 0.0;
    pbe[i] = Vpbe_ctor(alist[pbeparm->molid-1], pbeparm->nion,
      pbeparm->ionc, pbeparm->ionr, pbeparm->ionq, pbeparm->temp,
      pbeparm->gamma, pbeparm->pdie, pbeparm->sdie, sparm);

    /* Set up PDE object */
    pmgp[i] = Vpmgp_ctor(mgparm->dime[0], mgparm->dime[1],
      mgparm->dime[2], mgparm->nlev, mgparm->grid[0], mgparm->grid[1],
      mgparm->grid[2], pbeparm->nonlin);
    pmgp[i]->bcfl = pbeparm->bcfl;
    pmgp[i]->xcent = realCenter[0];
    pmgp[i]->ycent = realCenter[1];
    pmgp[i]->zcent = realCenter[2];
    if (pbeparm->bcfl == 4) {
        if (i == 0) {
            Vnm_tprint( 2, "Can't focus first calculation!\n");
            return 0;
        }
        pmg[i] = Vpmg_ctorFocus(pmgp[i], pbe[i], pmg[i-1],
          pbeparm->calcenergy);
    } else {
        if (i>0) Vpmg_dtor(&(pmg[i-1]));
        pmg[i] = Vpmg_ctor(pmgp[i], pbe[i]);
    }
    if (i>0) {
        Vpmgp_dtor(&(pmgp[i-1]));
        Vpbe_dtor(&(pbe[i-1]));
    }
    if (pbeparm->useDielMap) theDielXMap = dielXMap[pbeparm->dielMapID-1];
    else theDielXMap = VNULL;
    if (pbeparm->useDielMap) theDielYMap = dielYMap[pbeparm->dielMapID-1];
    else theDielYMap = VNULL;
    if (pbeparm->useDielMap) theDielZMap = dielZMap[pbeparm->dielMapID-1];
    else theDielZMap = VNULL;
    if (pbeparm->useKappaMap) theKappaMap = kappaMap[pbeparm->kappaMapID-1];
    else theKappaMap = VNULL;
    if (pbeparm->useChargeMap) theChargeMap = chargeMap[pbeparm->chargeMapID-1];
    else theChargeMap = VNULL;
    Vpmg_fillco(pmg[i], 
      pbeparm->srfm, pbeparm->swin, pbeparm->chgm,
      pbeparm->useDielMap, theDielXMap,
      pbeparm->useDielMap, theDielYMap,
      pbeparm->useDielMap, theDielZMap,
      pbeparm->useKappaMap, theKappaMap,
      pbeparm->useChargeMap, theChargeMap);

    /* Print a few derived parameters */
    Vnm_tprint(1, "  Debye length:  %g A\n", Vpbe_getDeblen(pbe[i]));

    /* Setup time statistics */
    Vnm_tstop(27, "Setup timer");

    /* Memory statistics */
    bytesTotal = Vmem_bytesTotal();
    highWater = Vmem_highWaterTotal();

#ifndef VAPBSQUIET
    Vnm_tprint( 1, "  Current memory usage:  %4.3f MB total, \
%4.3f MB high water\n", (double)(bytesTotal)/(1024.*1024.),
      (double)(highWater)/(1024.*1024.));
#endif


    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  killMG
//
// Purpose:  Kill MG structures
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void killMG(NOsh *nosh, Vpbe *pbe[NOSH_MAXCALC], 
  Vpmgp *pmgp[NOSH_MAXCALC], Vpmg *pmg[NOSH_MAXCALC]) {
    
    int i;

#ifndef VAPBSQUIET
    Vnm_tprint(1, "Destroying multigrid structures.\n");
#endif

    Vpbe_dtor(&(pbe[nosh->ncalc-1]));
    Vpmg_dtor(&(pmg[nosh->ncalc-1]));
    Vpmgp_dtor(&(pmgp[nosh->ncalc-1]));

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  solveMG
//
// Purpose:  Solve a PDE wth MG 
//
// Args:     type   MGparm::type
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int solveMG(NOsh *nosh, Vpmg *pmg, int type) {

    int nx, ny, nz, i;

   
    if (nosh != VNULL) {
        if (nosh->bogus) return 1;
    }

    Vnm_tstart(28, "Solver timer");


    if (type != 3) {
#ifndef VAPBSQUIET
        Vnm_tprint( 1,"  Solving PDE (see io.mc* for details)...\n");
#endif
        Vpmg_solve(pmg);
    } else {
        Vnm_tprint( 1,"  Skipping solve for mg-dummy run; zeroing \
solution array\n");
        nx = pmg->pmgp->nx;
        ny = pmg->pmgp->ny;
        nz = pmg->pmgp->nz;
        for (i=0; i<nx*ny*nz; i++) pmg->u[i] = 0.0;
    }
    Vnm_tstop(27, "Solver timer");

    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  setPartMG
//
// Purpose:  Set partition information for observables and I/I
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int setPartMG(NOsh *nosh, MGparm *mgparm, Vpmg *pmg) {

    int j;
    double partMin[3], partMax[3];

    if (nosh->bogus) return 1;

    if (mgparm->type == 2) {
        for (j=0; j<3; j++) {
            partMin[j] = mgparm->center[j] + mgparm->partDisjCenterShift[j]
              - 0.5*mgparm->partDisjLength[j];
            partMax[j] = mgparm->center[j] + mgparm->partDisjCenterShift[j]
              + 0.5*mgparm->partDisjLength[j];
        }
        Vnm_tprint(0, "Disj part lower corner = (%g, %g, %g)\n",
          partMin[0], partMin[1], partMin[2]);
        Vnm_tprint(0, "Disj part upper corner = (%g, %g, %g)\n",
          partMax[0], partMax[1], partMax[2]);
    } else {
        for (j=0; j<3; j++) {
            partMin[j] = mgparm->center[j] - 0.5*mgparm->glen[j];
            partMax[j] = mgparm->center[j] + 0.5*mgparm->glen[j];
        }
    }
    Vpmg_setPart(pmg, partMin, partMax, mgparm->partDisjOwnSide);


    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  energyMG
//
// Purpose:  Calculate and write out energies for MG calculation
//
// Args:     nosh       Holds input file information
//           pmg        Holds solution
//           icalc      Calculation index in nosh
//           totEnergy  set to total energy
//           qfEnergy   set to charge-phi energy
//           qmEnergy   set to mobile ion energy
//           dielEnergy set to dielectric energy
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int energyMG(NOsh *nosh, int icalc, Vpmg *pmg, 
  int *nenergy, double *totEnergy, double *qfEnergy, double *qmEnergy,
  double *dielEnergy) {

    Valist *alist;
    Vatom *atom;
    int i;
    double tenergy;
    MGparm *mgparm;
    PBEparm *pbeparm;
    int extEnergy;              /* When focusing, do we include energy 
                                 * contributions from outside the local 
                                 * partition? */

    mgparm = nosh->calc[icalc].mgparm;
    pbeparm = nosh->calc[icalc].pbeparm;

    if (mgparm->type == 2) extEnergy = 0;
    else extEnergy = 1;

    if (pbeparm->calcenergy == 1) {
        *nenergy = 1;
        /* Some processors don't count */
        if (nosh->bogus == 0) {
            *totEnergy = Vpmg_energy(pmg, extEnergy);
#ifndef VAPBSQUIET
            Vnm_tprint( 1, "  Total electrostatic energy = \
%1.12E kJ/mol\n", Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*totEnergy));
#endif
        } else *totEnergy = 0;
    } else if (pbeparm->calcenergy == 2) {
        *nenergy = 1;
        *totEnergy = Vpmg_energy(pmg, extEnergy);
        *qfEnergy = Vpmg_qfEnergy(pmg, extEnergy);
        *qmEnergy = Vpmg_qmEnergy(pmg, extEnergy);
        *dielEnergy = Vpmg_dielEnergy(pmg, extEnergy);
        Vnm_tprint( 1, "  Total electrostatic energy = %1.12E \
kJ/mol\n", Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*totEnergy));
        Vnm_tprint( 1, "  Fixed charge energy = %g kJ/mol\n",
           0.5*Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*qfEnergy));
        Vnm_tprint( 1, "  Mobile charge energy = %g kJ/mol\n",
           Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*qmEnergy));
        Vnm_tprint( 1, "  Dielectric energy = %g kJ/mol\n",
           Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*dielEnergy));
        Vnm_tprint( 1, "  Per-atom energies:\n");
        alist = pmg->pbe->alist;
        for (i=0; i<Valist_getNumberAtoms(alist); i++) {
            atom = Valist_getAtom(alist, i); 
            tenergy = Vpmg_qfAtomEnergy(pmg, atom);
            Vnm_tprint( 1, "      Atom %d:  %1.12E kJ/mol\n", i,
              0.5*Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*tenergy);
        }
    } else *nenergy = 0;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  forceMG
//
// Purpose:  Calculate and write out forces for MG calculation
//
// Args:     mem         Memory management
//           nosh        stores input file information
//           pbeparm     PBE parameters
//           pmg         Vpmg object for calculation
//           nforce      0 => no forces, 1 => net forces, >1 => number of
//                       forces (1 per atom)
//           atomForce   pointer to array of force objects
//           alist       molecules
//
// Returns:  1 if sucessful, 0 otherwise
// 
// Notes:    Sometimes (if nosh->bogus == 1) we just go through the motions,
//           but don't assign any forces
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int forceMG(Vmem *mem, NOsh *nosh, PBEparm *pbeparm, 
   Vpmg *pmg, int *nforce, AtomForce **atomForce, Valist *alist[NOSH_MAXMOL]) {

    int j, k;
    double qfForce[3], dbForce[3], ibForce[3], npForce[3];

    if (pbeparm->calcforce == 1) {
        *nforce = 1;
        *atomForce = (AtomForce *)Vmem_malloc(mem, 1, sizeof(AtomForce));
        /* Clear out force arrays */
        for (j=0; j<3; j++) {
            (*atomForce)[0].qfForce[j] = 0;
            (*atomForce)[0].ibForce[j] = 0;
            (*atomForce)[0].dbForce[j] = 0;
            (*atomForce)[0].npForce[j] = 0;
        }
        for (j=0;j<Valist_getNumberAtoms(alist[pbeparm->molid-1]);j++) { 
            if (nosh->bogus == 0) {
                Vpmg_qfForce(pmg, qfForce, j, pbeparm->chgm);
                Vpmg_ibForce(pmg, ibForce, j, pbeparm->srfm);
                Vpmg_dbnpForce(pmg, dbForce, npForce, j, pbeparm->srfm);
            } else {
                for (k=0; k<3; k++) {
                    qfForce[k] = 0; 
                    ibForce[k] = 0; 
                    dbForce[k] = 0; 
                    npForce[k] = 0; 
                }
            }
            for (k=0; k<3; k++) {
                (*atomForce)[0].qfForce[k] += qfForce[k];
                (*atomForce)[0].ibForce[k] += ibForce[k];
                (*atomForce)[0].dbForce[k] += dbForce[k];
                (*atomForce)[0].npForce[k] += npForce[k];
            }
        }
        Vnm_tprint( 1, "  Net fixed charge force on molecule %d\n",
          pbeparm->molid);
        Vnm_tprint( 1, "           = (%4.3e, %4.3e, %4.3e) kJ/(mol A)\n",
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].qfForce[0],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].qfForce[1],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].qfForce[2]);
        Vnm_tprint( 1, "  Net ionic boundary force on molecule %d\n",
          pbeparm->molid);
        Vnm_tprint( 1, "           = (%4.3e, %4.3e, %4.3e) kJ/(mol A)\n",
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].ibForce[0],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].ibForce[1],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].ibForce[2]);
        Vnm_tprint( 1, "  Net dielectric boundary force on \
molecule %d\n", pbeparm->molid);
        Vnm_tprint( 1, "           = (%4.3e, %4.3e, %4.3e) kJ/(mol A)\n",
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].dbForce[0],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].dbForce[1],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].dbForce[2]);
        Vnm_tprint( 1, "  Net apolar force on molecule %d\n",
          pbeparm->molid);
        Vnm_tprint( 1, "           = (%4.3e, %4.3e, %4.3e) kJ/(mol A)\n",
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].npForce[0],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].npForce[1],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].npForce[2]);
    } else if (pbeparm->calcforce == 2) {
        *nforce = Valist_getNumberAtoms(alist[pbeparm->molid-1]);
        *atomForce = (AtomForce *)Vmem_malloc(mem, *nforce,
          sizeof(AtomForce));
        for (j=0;j<Valist_getNumberAtoms(alist[pbeparm->molid-1]);j++) {
            if (nosh->bogus == 0) {
                Vpmg_qfForce(pmg, (*atomForce)[j].qfForce, j, pbeparm->chgm);
                Vpmg_ibForce(pmg, (*atomForce)[j].ibForce, j, pbeparm->srfm);
                Vpmg_dbnpForce(pmg, (*atomForce)[j].dbForce,
                  (*atomForce)[j].npForce, j, pbeparm->srfm);
            } else {
                for (k=0; k<3; k++) {
                    (*atomForce)[j].qfForce[k] = 0;
                    (*atomForce)[j].ibForce[k] = 0;
                    (*atomForce)[j].dbForce[k] = 0;
                    (*atomForce)[j].npForce[k] = 0;
                }
            }
            Vnm_tprint( 1, "  Total force on atom %d, molecule %d \
= (%4.3e, %4.3e, %4.3e) kJ/(mol A)\n", j, pbeparm->molid,
              Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(
                (*atomForce)[j].qfForce[0]+(*atomForce)[j].ibForce[0]+
                (*atomForce)[j].dbForce[0]+(*atomForce)[j].npForce[0]),
              Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(
                (*atomForce)[j].qfForce[1]+(*atomForce)[j].ibForce[1]+
                (*atomForce)[j].dbForce[1]+(*atomForce)[j].npForce[1]),
              Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(
                (*atomForce)[j].qfForce[2]+(*atomForce)[j].ibForce[2]+
                (*atomForce)[j].dbForce[2]+(*atomForce)[j].npForce[2]));
            Vnm_tprint( 1, "  Fixed charge force on atom %d, \
molecule %d = (%4.3e, %4.3e, %4.3e) kJ/mol/A\n", j, pbeparm->molid,
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].qfForce[0],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].qfForce[1],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].qfForce[2]);
            Vnm_tprint( 1, "  Ionic boundary force on atom %d, \
molecule %d = (%4.3e, %4.3e, %4.3e) kJ/mol/A\n", j, pbeparm->molid,
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].ibForce[0],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].ibForce[1],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].ibForce[2]);
            Vnm_tprint( 1, "  Dielectric boundary force on atom \
%d, molecule %d = (%4.3e, %4.3e, %4.3e) kJ/mol/A\n", j, pbeparm->molid,
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].dbForce[0],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].dbForce[1],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].dbForce[2]);
            Vnm_tprint( 1, "  Apolar force on atom %d, molecule %d \
= (%4.3e, %4.3e, %4.3e) kJ/mol/A\n", j, pbeparm->molid,
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].npForce[0],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].npForce[1],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].npForce[2]);
        }
    } else *nforce = 0;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  killEnergy
//
// Purpose:  Clear out energy structures 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void killEnergy() { 

#ifndef VAPBSQUIET
    Vnm_tprint(1, "No energy arrays to destroy.\n"); 
#endif

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  killForce
//
// Purpose:  Clear out force structures 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void killForce(Vmem *mem, NOsh *nosh, int nforce[NOSH_MAXCALC], 
  AtomForce *atomForce[NOSH_MAXCALC]) {

    int i;

#ifndef VAPBSQUIET
    Vnm_tprint(1, "Destroying force arrays.\n");
#endif

    for (i=0; i<nosh->ncalc; i++) {

        if (nforce[i] > 0) Vmem_free(mem, nforce[i], sizeof(AtomForce),
          (void **)&(atomForce[i]));
        
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  writematMG
//
// Purpose:  Write out matrix for MG calculation
//
// Returns:  1 if sucessful, 0 otherwise
//
// Notes:    currently ignores partition information when writing out acc
// 
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int writematMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg) {

    char writematstem[VMAX_ARGLEN];
    char outpath[VMAX_ARGLEN];
    char mxtype[3];
    int strlenmax;

    if (nosh->bogus) return 1;

#ifdef HAVE_MPI_H
    strlenmax = VMAX_ARGLEN-14;
    if (strlen(pbeparm->writematstem) > strlenmax) {
        Vnm_tprint(2, "  Matrix name (%s) too long (%d char max)!\n",
          pbeparm->writematstem, strlenmax);
        Vnm_tprint(2, "  Not writing matrix!\n");
        return 0;
    }
    sprintf(writematstem, "%s-PE%d", pbeparm->writematstem, rank);
#else
    strlenmax = VMAX_ARGLEN-1;
    if (strlen(pbeparm->writematstem) > strlenmax) {
        Vnm_tprint(2, "  Matrix name (%s) too long (%d char max)!\n",
          pbeparm->writematstem, strlenmax);
        Vnm_tprint(2, "  Not writing matrix!\n");
        return 0;
    }
    sprintf(writematstem, "%s", pbeparm->writematstem);
#endif
    
    if (pbeparm->writemat == 1) {
        strlenmax = VMAX_ARGLEN-5;
        if (strlen(pbeparm->writematstem) > strlenmax) {
            Vnm_tprint(2, "  Matrix name (%s) too long (%d char max)!\n",
              pbeparm->writematstem, strlenmax);
            Vnm_tprint(2, "  Not writing matrix!\n");
            return 0;
        }
        sprintf(outpath, "%s.%s", writematstem, "mat");
        mxtype[0] = 'R';
        mxtype[1] = 'S';
        mxtype[2] = 'A';
        /* Poisson operator only */
        if (pbeparm->writematflag == 0) {
            Vnm_tprint( 1, "  Writing Poisson operator matrix \
to %s...\n", outpath);

         /* Linearization of Poisson-Boltzmann operator around solution */
         } else if (pbeparm->writematflag == 1) {
            Vnm_tprint( 1, "  Writing linearization of full \
Poisson-Boltzmann operator matrix to %s...\n", outpath);

         } else {
             Vnm_tprint( 2, "  Bogus matrix specification\
(%d)!\n", pbeparm->writematflag);
             return 0;
         }

         Vnm_tprint(0, "  Printing operator...\n");
         Vpmg_printColComp(pmg, outpath, outpath, mxtype, 
           pbeparm->writematflag);

    }

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  writedataMG
//
// Purpose:  Write out data from  MG calculation
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int writedataMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg) {

    char writestem[VMAX_ARGLEN];
    char outpath[VMAX_ARGLEN];
    char title[72];
    int i, nx, ny, nz;
    double hx, hy, hzed, xcent, ycent, zcent, xmin, ymin, zmin;
    Vgrid *grid; 

    if (nosh->bogus) return 1;

  
    for (i=0; i<pbeparm->numwrite; i++) { 

        nx = pmg->pmgp->nx;
        ny = pmg->pmgp->ny;
        nz = pmg->pmgp->nz;
        hx = pmg->pmgp->hx;
        hy = pmg->pmgp->hy;
        hzed = pmg->pmgp->hzed;

        switch (pbeparm->writetype[i]) {

            case VDT_CHARGE:
 
                Vnm_tprint(1, "  Writing charge distribution to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_CHARGE, 0.0);
                sprintf(title, "CHARGE DISTRIBUTION (e)");
                break;

            case VDT_POT:
 
                Vnm_tprint(1, "  Writing potential to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_POT, 0.0);
                sprintf(title, "POTENTIAL (kT/e)");
                break;

            case VDT_SMOL:

                Vnm_tprint(1, "  Writing molecular accessibility to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_SMOL, pbeparm->srad);
                sprintf(title, 
                  "SOLVENT ACCESSIBILITY -- MOLECULAR (%4.3f PROBE)", 
                  pbeparm->srad);
                break;

            case VDT_SSPL:

                Vnm_tprint(1, "  Writing spline-based accessibility to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_SSPL, pbeparm->swin);
                sprintf(title, 
                  "SOLVENT ACCESSIBILITY -- SPLINE (%4.3f WINDOW)",
                  pbeparm->swin);
                break;

            case VDT_VDW:

                Vnm_tprint(1, "  Writing van der Waals accessibility to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_VDW, 0.0);
                sprintf(title, "SOLVENT ACCESSIBILITY -- VAN DER WAALS");
                break;

            case VDT_IVDW:

                Vnm_tprint(1, "  Writing ion accessibility to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_IVDW, 
                  pmg->pbe->maxIonRadius);
                sprintf(title, 
                  "ION ACCESSIBILITY -- SPLINE (%4.3f RADIUS)",
                  pmg->pbe->maxIonRadius);
                break;

            case VDT_LAP:

                Vnm_tprint(1, "  Writing potential Laplacian to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_LAP, 0.0);
                sprintf(title, 
                  "POTENTIAL LAPLACIAN (kT/e/A^2)");
                break;

            case VDT_EDENS:

                Vnm_tprint(1, "  Writing energy density to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_EDENS, 0.0);
                sprintf(title, "ENERGY DENSITY (kT/e/A)^2");
                break;

            case VDT_NDENS:

                Vnm_tprint(1, "  Writing number density to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_NDENS, 0.0);
                sprintf(title, 
                  "ION NUMBER DENSITY (M)");
                break;

            case VDT_QDENS:

                Vnm_tprint(1, "  Writing charge density to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_QDENS, 0.0);
                sprintf(title, 
                  "ION CHARGE DENSITY (e_c * M)");
                break;

            case VDT_DIELX:

                Vnm_tprint(1, "  Writing x-shifted dielectric map to ");
                xcent = pmg->pmgp->xcent + 0.5*hx;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_DIELX, 0.0);
                sprintf(title,
                  "X-SHIFTED DIELECTRIC MAP");
                break;

            case VDT_DIELY:

                Vnm_tprint(1, "  Writing y-shifted dielectric map to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent + 0.5*hy;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_DIELY, 0.0);
                sprintf(title,
                  "Y-SHIFTED DIELECTRIC MAP");
                break;

            case VDT_DIELZ:

                Vnm_tprint(1, "  Writing z-shifted dielectric map to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent + 0.5*hzed;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_DIELZ, 0.0);
                sprintf(title,
                  "Z-SHIFTED DIELECTRIC MAP");
                break;

            case VDT_KAPPA:

                Vnm_tprint(1, "  Writing kappa map to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_KAPPA, 0.0);
                sprintf(title,
                  "KAPPA MAP");
                break;

            default:

                Vnm_tprint(2, "Invalid data type for writing!\n");
                return 0;
                break;
        }


#ifdef HAVE_MPI_H
        sprintf(writestem, "%s-PE%d", pbeparm->writestem[i], rank);
#else
        sprintf(writestem, "%s", pbeparm->writestem[i]);
#endif

        switch (pbeparm->writefmt[i]) {

            case VDF_DX:
                sprintf(outpath, "%s.%s", writestem, "dx");
                Vnm_tprint(1, "%s\n", outpath);
                grid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin,
                  pmg->rwork);
                Vgrid_writeDX(grid, "FILE", "ASC", VNULL, outpath, title,
                  pmg->pvec);
                Vgrid_dtor(&grid);
                break;

            case VDF_AVS:
                sprintf(outpath, "%s.%s", writestem, "ucd");
                Vnm_tprint(1, "%s\n", outpath);
                Vnm_tprint(2, "Sorry, AVS format isn't supported for \
uniform meshes yet!\n");
                break;

            case VDF_UHBD:
                sprintf(outpath, "%s.%s", writestem, "grd");
                Vnm_tprint(1, "%s\n", outpath);
                grid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin,
                  pmg->rwork);
                Vgrid_writeUHBD(grid, "FILE", "ASC", VNULL, outpath, title,
                  pmg->pvec);
                Vgrid_dtor(&grid);
                break;

            default:
                Vnm_tprint(2, "Bogus data format (%d)!\n", 
                  pbeparm->writefmt[i]);
                break;
        }
                
    }

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  printEnergy
//
// Purpose:  Execute a PRINT ENERGY statement
//
// Args:     i     Index of energy statement to print
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int printEnergy(Vcom *com, NOsh *nosh, double totEnergy[NOSH_MAXCALC], 
  int i) {

    int j, calcid;
    double ltenergy, gtenergy, scalar;

    Vnm_tprint( 1, "print energy %d ", nosh->printcalc[i][0]);
    for (j=1; j<nosh->printnarg[i]; j++) {
        if (nosh->printop[i][j-1] == 0)
          Vnm_tprint(1, "+ ");
        else if (nosh->printop[i][j-1] == 1)
          Vnm_tprint(1, "- ");
        else {
            Vnm_tprint( 2, "Undefined PRINT operation!\n");
            return 0;
        }
        Vnm_tprint(1, "%d ", nosh->printcalc[i][j]);
    }
    Vnm_tprint(1, "end\n");
    calcid = nosh->elec2calc[nosh->printcalc[i][0]-1];
    if (nosh->calc[calcid].pbeparm->calcenergy != 0) {
        ltenergy = Vunit_kb * (1e-3) * Vunit_Na *
          nosh->calc[calcid].pbeparm->temp * totEnergy[calcid];
    } else {
        Vnm_tprint( 2, "  Didn't calculate energy in Calculation \
#%d\n", calcid+1);
        return 0;
    }
    for (j=1; j<nosh->printnarg[i]; j++) {
        calcid = nosh->elec2calc[nosh->printcalc[i][j]-1];
        /* Add or subtract? */
        if (nosh->printop[i][j-1] == 0) scalar = 1.0;
        else if (nosh->printop[i][j-1] == 1) scalar = -1.0;
        /* Accumulate */
        ltenergy += (scalar * Vunit_kb * (1e-3) * Vunit_Na *
          nosh->calc[calcid].pbeparm->temp * totEnergy[calcid]);
    }

    Vnm_tprint( 1, "  Local net energy (PE %d) = %1.12E kJ/mol\n", 
      Vcom_rank(com), ltenergy);
    Vnm_tprint( 0, "printEnergy:  Performing global reduction (sum)\n");
    Vcom_reduce(com, &ltenergy, &gtenergy, 1, 2, 0);
    Vnm_tprint( 1, "  Global net energy = %1.12E kJ/mol\n", gtenergy);

    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  printForce
//
// Purpose:  Execute a PRINT ENERGY statement
//
// Args:     i     Index of energy statement to print
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int printForce(Vcom *com, NOsh *nosh, int nforce[NOSH_MAXCALC], 
  AtomForce *atomForce[NOSH_MAXCALC], int i) {

    int ipr, ifr, ivc, calcid, refnforce, refcalcforce;
    double temp, scalar;
    AtomForce *lforce, *gforce, *aforce;

    Vnm_tprint( 1, "print force %d ", nosh->printcalc[i][0]);
    for (ipr=1; ipr<nosh->printnarg[i]; ipr++) {
        if (nosh->printop[i][ipr-1] == 0)
          Vnm_tprint(1, "+ ");
        else if (nosh->printop[i][ipr-1] == 1)
          Vnm_tprint(1, "- ");
        else {
            Vnm_tprint( 2, "Undefined PRINT operation!\n");
            return 0;
        }
        Vnm_tprint(1, "%d ", nosh->printcalc[i][ipr]);
    }
    Vnm_tprint(1, "end\n");

    /* First, go through and make sure we did the same type of force
     * evaluation in each of the requested calculations */
    calcid = nosh->elec2calc[nosh->printcalc[i][0]-1];
    refnforce = nforce[calcid];
    refcalcforce = nosh->calc[calcid].pbeparm->calcforce;
    if (refcalcforce == 0) {
        Vnm_tprint( 2, "  Didn't calculate force in calculation \
#%d\n", calcid+1);
        return 0;
    }
    for (ipr=1; ipr<nosh->printnarg[i]; ipr++) {
        calcid = nosh->elec2calc[nosh->printcalc[i][ipr]-1];
        if (nosh->calc[calcid].pbeparm->calcforce != refcalcforce) {
            Vnm_tprint(2, "  Inconsistent calcforce declarations in \
calculations %d and %d\n", nosh->elec2calc[nosh->printcalc[i][0]-1]+1,
calcid+1);
            return 0;
        }
        if (nforce[calcid] != refnforce) {
            Vnm_tprint(2, "  Inconsistent number of forces evaluated in \
calculations %d and %d\n", nosh->elec2calc[nosh->printcalc[i][0]-1]+1,
calcid+1);
            return 0;
        }
    }

    /* Now, allocate space to accumulate the forces */
    lforce = (AtomForce *)Vmem_malloc(VNULL, refnforce, sizeof(AtomForce));
    gforce = (AtomForce *)Vmem_malloc(VNULL, refnforce, sizeof(AtomForce));

    /* Now, accumulate the forces */
    calcid = nosh->elec2calc[nosh->printcalc[i][0]-1];
    aforce = atomForce[calcid];
    temp = nosh->calc[calcid].pbeparm->temp;

    /* Load up the first calculation */
    if (refcalcforce == 1) {
        /* Set to total force */
        for (ivc=0; ivc<3; ivc++) {
            lforce[0].qfForce[ivc] = 
              Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].qfForce[ivc];
            lforce[0].ibForce[ivc] = 
              Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].ibForce[ivc];
            lforce[0].dbForce[ivc] = 
              Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].dbForce[ivc];
            lforce[0].npForce[ivc] = 
              Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].npForce[ivc];
        }
    } else if (refcalcforce == 2) { 
        for (ifr=0; ifr<refnforce; ifr++) {
            for (ivc=0; ivc<3; ivc++) {
                lforce[ifr].qfForce[ivc] = 
              Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].qfForce[ivc];
                lforce[ifr].ibForce[ivc] = 
              Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].ibForce[ivc];
                lforce[ifr].dbForce[ivc] = 
              Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].dbForce[ivc];
                lforce[ifr].npForce[ivc] = 
              Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].npForce[ivc];
            }
        }
    }

    /* Load up the rest of the calculations */
    for (ipr=1; ipr<nosh->printnarg[i]; ipr++) {
        calcid = nosh->elec2calc[nosh->printcalc[i][ipr]-1];
        temp = nosh->calc[calcid].pbeparm->temp;
        aforce = atomForce[calcid];
        /* Get operation */
        if (nosh->printop[i][ipr-1] == 0) scalar = +1.0;
        else if (nosh->printop[i][ipr-1] == 1) scalar = -1.0;
        else scalar = 0.0;
        /* Accumulate */
        if (refcalcforce == 1) {
            /* Set to total force */
            for (ivc=0; ivc<3; ivc++) {
                lforce[0].qfForce[ivc] += 
                 (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].qfForce[ivc]);
                lforce[0].ibForce[ivc] += 
                 (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].ibForce[ivc]);
                lforce[0].dbForce[ivc] += 
                 (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].dbForce[ivc]);
                lforce[0].npForce[ivc] +=
                 (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].npForce[ivc]);
            }
        } else if (refcalcforce == 2) {
            for (ifr=0; ifr<refnforce; ifr++) {
                for (ivc=0; ivc<3; ivc++) {
                    lforce[ifr].qfForce[ivc] += 
               (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].qfForce[ivc]);
                    lforce[ifr].ibForce[ivc] += 
               (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].ibForce[ivc]);
                    lforce[ifr].dbForce[ivc] += 
               (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].dbForce[ivc]);
                    lforce[ifr].npForce[ivc] += 
               (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].npForce[ivc]);
                }
            }
        }
    }

    Vnm_tprint( 0, "printEnergy:  Performing VERY INEFFICIENT global reduction (sum)\n");
    for (ifr=0; ifr<refnforce; ifr++) {
        Vcom_reduce(com, lforce[ifr].qfForce, gforce[ifr].qfForce, 3, 2, 0);
        Vcom_reduce(com, lforce[ifr].ibForce, gforce[ifr].ibForce, 3, 2, 0);
        Vcom_reduce(com, lforce[ifr].dbForce, gforce[ifr].dbForce, 3, 2, 0);
        Vcom_reduce(com, lforce[ifr].npForce, gforce[ifr].npForce, 3, 2, 0);
    }
   
#if 1
    if (refcalcforce == 1) {
        Vnm_tprint( 1, "  Local net fixed charge force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", lforce[0].qfForce[0],
lforce[0].qfForce[1], lforce[0].qfForce[2]);
        Vnm_tprint( 1, "  Local net ionic boundary force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", lforce[0].ibForce[0],
lforce[0].ibForce[1], lforce[0].ibForce[2]);
        Vnm_tprint( 1, "  Local net dielectric boundary force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", lforce[0].dbForce[0],
lforce[0].dbForce[1], lforce[0].dbForce[2]);
        Vnm_tprint( 1, "  Local net apolar boundary force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", lforce[0].npForce[0],
lforce[0].npForce[1], lforce[0].npForce[2]);
    } else if (refcalcforce == 2) {
        for (ifr=0; ifr<refnforce; ifr++) {
            Vnm_tprint( 1, "  Local fixed charge force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, lforce[ifr].qfForce[0],
lforce[ifr].qfForce[1], lforce[ifr].qfForce[2]);
        Vnm_tprint( 1, "  Local ionic boundary force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, lforce[ifr].ibForce[0],
lforce[ifr].ibForce[1], lforce[ifr].ibForce[2]);
        Vnm_tprint( 1, "  Local dielectric boundary force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, lforce[ifr].dbForce[0],
lforce[ifr].dbForce[1], lforce[ifr].dbForce[2]);
        Vnm_tprint( 1, "  Local apolar boundary force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, lforce[ifr].npForce[0],
lforce[ifr].npForce[1], lforce[ifr].npForce[2]);
        }
    }
#endif
 
    if (refcalcforce == 1) {
        Vnm_tprint( 1, "  Global net fixed charge force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", gforce[0].qfForce[0], 
gforce[0].qfForce[1], gforce[0].qfForce[2]);
        Vnm_tprint( 1, "  Global net ionic boundary force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", gforce[0].ibForce[0], 
gforce[0].ibForce[1], gforce[0].ibForce[2]);
        Vnm_tprint( 1, "  Global net dielectric boundary force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", gforce[0].dbForce[0], 
gforce[0].dbForce[1], gforce[0].dbForce[2]);
        Vnm_tprint( 1, "  Global net apolar boundary force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", gforce[0].npForce[0], 
gforce[0].npForce[1], gforce[0].npForce[2]);
    } else if (refcalcforce == 2) {
        for (ifr=0; ifr<refnforce; ifr++) {
            Vnm_tprint( 1, "  Global fixed charge force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, gforce[ifr].qfForce[0], 
gforce[ifr].qfForce[1], gforce[ifr].qfForce[2]);
        Vnm_tprint( 1, "  Global ionic boundary force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, gforce[ifr].ibForce[0],
gforce[ifr].ibForce[1], gforce[ifr].ibForce[2]);
        Vnm_tprint( 1, "  Global dielectric boundary force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, gforce[ifr].dbForce[0],
gforce[ifr].dbForce[1], gforce[ifr].dbForce[2]);
        Vnm_tprint( 1, "  Global apolar boundary force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, gforce[ifr].npForce[0],
gforce[ifr].npForce[1], gforce[ifr].npForce[2]);
        Vnm_tprint( 1, "  Global total force (atom %d) = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, 
(gforce[ifr].npForce[0] + gforce[ifr].dbForce[0] + gforce[ifr].ibForce[0] +
gforce[ifr].qfForce[0]),
(gforce[ifr].npForce[1] + gforce[ifr].dbForce[1] + gforce[ifr].ibForce[1] +
gforce[ifr].qfForce[1]),
(gforce[ifr].npForce[2] + gforce[ifr].dbForce[2] + gforce[ifr].ibForce[2] +
gforce[ifr].qfForce[2]));
        }
    }

    Vmem_free(VNULL, refnforce, sizeof(AtomForce), (void **)(&lforce));
    Vmem_free(VNULL, refnforce, sizeof(AtomForce), (void **)(&gforce));

    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  npenergyMG
//
// Purpose:  Calculate and write out energies for MG calculation
//
// Args:     nosh       Holds input file information
//           pmg        Holds solution
//           icalc      Calculation index in nosh
//           npEnergy   set to apolar energy
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Robert Konecny (based on energyMG)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int npenergyMG(NOsh *nosh, int icalc, Vpmg *pmg,
  int *nenergy, double *npEnergy) {

    MGparm *mgparm;
    PBEparm *pbeparm;
    int extEnergy;              /* When focusing, do we include energy
                                 * contributions from outside the local
                                 * partition? */

    mgparm = nosh->calc[icalc].mgparm;
    pbeparm = nosh->calc[icalc].pbeparm;

    if (mgparm->type == 2) extEnergy = 0;
    else extEnergy = 1;

    if (pbeparm->calcenergy > 0) {
        *nenergy = 1;
        /* Some processors don't count */
        if (nosh->bogus == 0) {
        *npEnergy = Vpmg_npEnergy(pmg, extEnergy);
        } else *npEnergy = 0;

    } else *nenergy = 0;

    return 1;
}
