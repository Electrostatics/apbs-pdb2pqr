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
    
    int i;

    Vnm_tprint( 1, "main:  Got PQR paths for %d molecules\n", nosh->nmol);
    if (nosh->nmol <= 0) {
       Vnm_tprint(2, "main:  You didn't specify any molecules (correctly)!\n");
       Vnm_tprint(2, "main:  Bailing out!\n");
       return 0;
    }

    for (i=0; i<nosh->nmol; i++) {
        Vnm_tprint( 1, "main:  Reading atom data from %s:\n",
          nosh->molpath[i]);
        alist[i] = Valist_ctor();
        if (Valist_readPQR(alist[i], "FILE", "ASC", VNULL,
          nosh->molpath[i]) != 1) {
            Vnm_tprint( 2, "main:  Fatal error while reading from %s\n",
              nosh->molpath[i]);
            return 0;
        } else {
            Vnm_tprint( 1, "main:    %d atoms\n",
              Valist_getNumberAtoms(alist[i]));
            Vnm_tprint( 1, "main:    Centered at (%4.3e, %4.3e, %4.3e)\n",
              alist[i]->center[0], alist[i]->center[1], alist[i]->center[2]);
            Vnm_tprint( 1, "main:    Net charge %3.2e e\n",
              alist[i]->charge);        
        }
    }

    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  loadDielMaps
// 
// Purpose:  Load dielectric maps from files
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
      Vnm_tprint( 1, "main:  Got paths for %d dielectric map sets\n", 
        nosh->ndiel);
    else return 1;

    for (i=0; i<nosh->ndiel; i++) {
        Vnm_tprint( 1, "main:  Reading x-shifted dielectric map data from \
%s:\n", nosh->dielXpath[i]);
        dielXMap[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        if (nosh->dielfmt[i] == 0) {
            if (Vgrid_readDX(dielXMap[i], "FILE", "ASC", VNULL, 
              nosh->dielXpath[i]) != 1) {
                Vnm_tprint( 2, "main:  Fatal error while reading from %s\n",
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
            Vnm_tprint(1, "main:    %d x %d x %d grid\n", nx, ny, nz);
            Vnm_tprint(1, "main:    (%g, %g, %g) A spacings\n", hx, hy, hzed);
            Vnm_tprint(1, "main:    (%g, %g, %g) A lower corner\n", 
              xmin, ymin, zmin);
            sum = 0;
            for (ii=0; ii<(nx*ny*nz); ii++)
              sum += (dielXMap[i]->data[ii]);
            sum = sum*hx*hy*hzed;
            Vnm_print(2, "main:    Volume integral = %3.2e A^3\n", 
              sum);
        } else {
            Vnm_tprint( 2, "main:  INVALID FORMAT!\n");
            return 0;
        }
        Vnm_tprint( 1, "main:  Reading y-shifted dielectric map data from \
%s:\n", nosh->dielYpath[i]);
        dielYMap[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        if (nosh->dielfmt[i] == 0) {
            if (Vgrid_readDX(dielYMap[i], "FILE", "ASC", VNULL, nosh->dielYpath[i])
              != 1) {
                Vnm_tprint( 2, "main:  Fatal error while reading from %s\n",
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
            Vnm_tprint(1, "main:    %d x %d x %d grid\n", nx, ny, nz);
            Vnm_tprint(1, "main:    (%g, %g, %g) A spacings\n", hx, hy, hzed);
            Vnm_tprint(1, "main:    (%g, %g, %g) A lower corner\n",
              xmin, ymin, zmin);
            sum = 0;
            for (ii=0; ii<(nx*ny*nz); ii++)
              sum += (dielYMap[i]->data[ii]);
            sum = sum*hx*hy*hzed;
            Vnm_print(2, "main:    Volume integral = %3.2e A^3\n",
              sum);
        } else {
            Vnm_tprint( 2, "main:  INVALID FORMAT!\n");
            return 0;
        }
        Vnm_tprint( 1, "main:  Reading z-shifted dielectric map data from \
%s:\n", nosh->dielZpath[i]);
        dielZMap[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        if (nosh->dielfmt[i] == 0) {
            if (Vgrid_readDX(dielZMap[i], "FILE", "ASC", VNULL, nosh->dielZpath[i])
              != 1) {
                Vnm_tprint( 2, "main:  Fatal error while reading from %s\n",
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
            Vnm_tprint(1, "main:    %d x %d x %d grid\n",
              nx, ny, nz);
            Vnm_tprint(1, "main:    (%g, %g, %g) A spacings\n",
              hx, hy, hzed);
            Vnm_tprint(1, "main:    (%g, %g, %g) A lower corner\n",
              xmin, ymin, zmin);
            sum = 0;
            for (ii=0; ii<(nx*ny*nz); ii++)
              sum += (dielZMap[i]->data[ii]);
            sum = sum*hx*hy*hzed;
            Vnm_print(2, "main:    Volume integral = %3.2e A^3\n",
              sum);
        } else {
            Vnm_tprint( 2, "main:  INVALID FORMAT!\n");
            return 0;
        }
    }

    return 1;

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
      Vnm_tprint( 1, "main:  Got paths for %d kappa maps\n", nosh->nkappa);
    else return 1;

    for (i=0; i<nosh->nkappa; i++) {
        Vnm_tprint( 1, "main:  Reading kappa map data from %s:\n",
          nosh->kappapath[i]);
        map[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        if (nosh->kappafmt[i] == 0) {
            if (Vgrid_readDX(map[i], "FILE", "ASC", VNULL, nosh->kappapath[i]) 
              != 1) {
                Vnm_tprint( 2, "main:  Fatal error while reading from %s\n",
                  nosh->kappapath[i]);
                return 0;
            }
            Vnm_tprint(1, "main:    %d x %d x %d grid\n", 
              map[i]->nx, map[i]->ny, map[i]->nz);
            Vnm_tprint(1, "main:    (%g, %g, %g) A spacings\n", 
              map[i]->hx, map[i]->hy, map[i]->hzed);
            Vnm_tprint(1, "main:    (%g, %g, %g) A lower corner\n", 
              map[i]->xmin, map[i]->ymin, map[i]->zmin);
            sum = 0;
            for (ii=0; ii<(map[i]->nx*map[i]->ny*map[i]->nz); ii++)
              sum += (map[i]->data[ii]);
            sum = sum*map[i]->hx*map[i]->hy*map[i]->hzed;
            Vnm_print(2, "main:    Volume integral = %3.2e A^3\n", sum);

        } else {
            Vnm_tprint( 2, "main:  INVALID FORMAT!\n");
            return 0;
        }
    }

    return 1;

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
      Vnm_tprint( 1, "main:  Got paths for %d charge maps\n", nosh->ncharge);
    else return 1;

    for (i=0; i<nosh->ncharge; i++) {
        Vnm_tprint( 1, "main:  Reading charge map data from %s:\n",
          nosh->chargepath[i]);
        map[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        if (nosh->chargefmt[i] == 0) {
            if (Vgrid_readDX(map[i], "FILE", "ASC", VNULL, nosh->chargepath[i])
              != 1) {
                Vnm_tprint( 2, "main:  Fatal error while reading from %s\n",
                  nosh->chargepath[i]);
                return 0;
            }
            Vnm_tprint(1, "main:    %d x %d x %d grid\n", 
              map[i]->nx, map[i]->ny, map[i]->nz);
            Vnm_tprint(1, "main:    (%g, %g, %g) A spacings\n", 
              map[i]->hx, map[i]->hy, map[i]->hzed);
            Vnm_tprint(1, "main:    (%g, %g, %g) A lower corner\n", 
              map[i]->xmin, map[i]->ymin, map[i]->zmin);
            sum = 0;
            for (ii=0; ii<(map[i]->nx*map[i]->ny*map[i]->nz); ii++) 
              sum += (map[i]->data[ii]);
            sum = sum*map[i]->hx*map[i]->hy*map[i]->hzed;
            Vnm_print(2, "main:    Charge map integral = %3.2e e\n", sum);
        } else {
            Vnm_tprint( 2, "main:  INVALID FORMAT!\n");
            return 0;
        }
    }

    return 1;

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

    Vnm_tprint( 1, "main:    Molecule ID: %d\n", pbeparm->molid);
    if (pbeparm->nonlin) Vnm_tprint( 1, "main:    Nonlinear PBE\n");
    else Vnm_tprint( 1, "main:    Linearized PBE\n");
    if (pbeparm->bcfl == 0) {
        Vnm_tprint( 1, "main:    Zero boundary conditions\n");
    } else if (pbeparm->bcfl == 1) {
        Vnm_tprint( 1, "main:    Single Debye-Huckel sphere boundary \
conditions\n");
    } else if (pbeparm->bcfl == 2) {
        Vnm_tprint( 1, "main:    Multiple Debye-Huckel sphere boundary \
conditions\n");
    } else if (pbeparm->bcfl == 4) {
        Vnm_tprint( 1, "main:    Boundary conditions from focusing\n");
    }
    Vnm_tprint( 1, "main:    %d ion species (%4.3f M ionic strength):\n",
      pbeparm->nion, ionstr);
    for (i=0; i<pbeparm->nion; i++) {
        Vnm_tprint( 1, "main:      %4.3f A-radius, %4.3f e-charge, \
%4.3f M concentration\n", 
          pbeparm->ionr[i], pbeparm->ionq[i], pbeparm->ionc[i]);            
    }
    Vnm_tprint( 1, "main:    Solute dielectric: %4.3f\n", pbeparm->pdie);
    Vnm_tprint( 1, "main:    Solvent dielectric: %4.3f\n", pbeparm->sdie);
    if (pbeparm->srfm == 0) {
        Vnm_tprint( 1, "main:    Using \"molecular\" surface \
definition; no smoothing\n");
        Vnm_tprint( 1, "main:    Solvent probe radius: %4.3f A\n",
          pbeparm->srad);
    } else if (pbeparm->srfm == 1) {
        Vnm_tprint( 1, "main:    Using \"molecular\" surface definition;\
 harmonic average smoothing\n");
        Vnm_tprint( 1, "main:    Solvent probe radius: %4.3f A\n",
          pbeparm->srad);
    } else if (pbeparm->srfm == 2) {
        Vnm_tprint( 1, "main:    Using spline-based surface definition;\
 window = %4.3f\n", pbeparm->swin);
    }
    Vnm_tprint( 1, "main:    Temperature:  %4.3f K\n", pbeparm->temp);
    Vnm_tprint( 1, "main:    Surface tension:  %4.3f kJ/mol/A^2\n",
      pbeparm->gamma);
    if (pbeparm->calcenergy == 1) Vnm_tprint( 1, "main:    Electrostatic \
energies will be calculated\n");
    if (pbeparm->calcforce == 1) Vnm_tprint( 1, "main:    Net solvent \
forces will be calculated \n");
    if (pbeparm->calcforce == 2) Vnm_tprint( 1, "main:    All-atom \
solvent forces will be calculated\n");
    for (i=0; i<pbeparm->numwrite; i++) {
        switch (pbeparm->writetype[i]) {
            case VDT_CHARGE:
                Vnm_print(1, "main:    Charge distribution to be written to ");
                break;
            case VDT_POT:
                Vnm_print(1, "main:    Potential to be written to ");
                break;
            case VDT_SMOL:
                Vnm_print(1, "main:    Molecular solvent accessibility \
to be written to ");
                break;
            case VDT_SSPL:
                Vnm_print(1, "main:    Spline-based solvent accessibility \
to be written to ");
                break;
            case VDT_VDW:
                Vnm_print(1, "main:    van der Waals solvent accessibility \
to be written to ");
                break;
            case VDT_IVDW:
                Vnm_print(1, "main:    Ion accessibility to be written to ");
                break;
            case VDT_LAP:
                Vnm_print(1, "main:    Potential Laplacian to be written to ");
                break;
            case VDT_EDENS:
                Vnm_print(1, "main:    Energy density to be written to ");
                break;
            case VDT_NDENS:
                Vnm_print(1, "main:    Ion number density to be written to ");
                break;
            case VDT_QDENS:
                Vnm_print(1, "main:    Ion charge density to be written to ");
                break;
            case VDT_DIELX:
                Vnm_print(1, "main:    X-shifted dielectric map to be written \
to ");
                break;
            case VDT_DIELY:
                Vnm_print(1, "main:    Y-shifted dielectric map to be written \
to ");
                break;
            case VDT_DIELZ:
                Vnm_print(1, "main:    Z-shifted dielectric map to be written \
to ");
                break;
            case VDT_KAPPA:
                Vnm_print(1, "main:    Kappa map to be written to ");
                break;
            default: 
                Vnm_print(2, "main:    Invalid data type for writing!\n");
                break;
        }
        switch (pbeparm->writefmt[i]) {
            case VDF_DX:
                Vnm_print(1, "%s.%s\n", pbeparm->writestem[i], "dx");
                break;
            case VDF_UHBD:
                Vnm_print(1, "%s.%s\n", pbeparm->writestem[i], "grd");
                break;
            case VDF_AVS:
                Vnm_print(1, "%s.%s\n", pbeparm->writestem[i], "ucd");
                break;
            default: 
                Vnm_print(2, "main:    Invalid format for writing!\n");
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
        Vnm_tprint( 1, "main:    Partition overlap fraction = %g\n", 
          mgparm->ofrac);
        Vnm_tprint( 1, "main:    Processor array = %d x %d x %d\n", 
          mgparm->pdime[0], mgparm->pdime[1], mgparm->pdime[2]);
    }
    Vnm_tprint( 1, "main:    Grid dimensions: %d x %d x %d\n",
      mgparm->dime[0], mgparm->dime[1], mgparm->dime[2]);
    Vnm_tprint( 1, "main:    Grid spacings: %4.3f x %4.3f x %4.3f\n",
      mgparm->grid[0], mgparm->grid[1], mgparm->grid[2]);
    Vnm_tprint( 1, "main:    Grid lengths: %4.3f x %4.3f x %4.3f\n",
      mgparm->glen[0], mgparm->glen[1], mgparm->glen[2]);
    Vnm_tprint( 1, "main:    Grid center: (%4.3f, %4.3f, %4.3f)\n",
      realCenter[0], realCenter[1], realCenter[2]);
    Vnm_tprint( 1, "main:    Multigrid levels: %d\n", mgparm->nlev);

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
    
    int j, bytesTotal, highWater;
    double sparm, iparm;
    Vgrid *theDielXMap, *theDielYMap, *theDielZMap, *theKappaMap, *theChargeMap;

    Vnm_tstart(27, "Setup timer");

    /* Fix mesh center for "GCENT MOL #" types of declarations */
    if (mgparm->cmeth == 1) {
        for (j=0; j<3; j++) 
          mgparm->center[j] = (alist[mgparm->centmol-1])->center[j];
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
      pbeparm->pdie, pbeparm->sdie, sparm);

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
            Vnm_tprint( 2, "main:  Can't focus first calculation!\n");
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
      pbeparm->srfm, pbeparm->swin,
      pbeparm->useDielMap, theDielXMap,
      pbeparm->useDielMap, theDielYMap,
      pbeparm->useDielMap, theDielZMap,
      pbeparm->useKappaMap, theKappaMap,
      pbeparm->useChargeMap, theChargeMap);

    /* Setup time statistics */
    Vnm_tstop(27, "Setup timer");

    /* Memory statistics */
    bytesTotal = Vmem_bytesTotal();
    highWater = Vmem_highWaterTotal();
    Vnm_tprint( 1, "main:    Current memory usage:  %4.3f MB total, \
%4.3f MB high water\n", (double)(bytesTotal)/(1024.*1024.),
      (double)(highWater)/(1024.*1024.));


    return 1;

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
VPUBLIC int solveMG(Vpmg *pmg, int type) {

    int nx, ny, nz, i;

    Vnm_tstart(28, "Solver timer");
    if (type != 3) {
        Vnm_tprint( 1,"main:    Solving PDE (see io.mc* for details)...\n");
        Vpmg_solve(pmg);
    } else {
        Vnm_tprint( 1,"main:    Skipping solve for mg-dummy run; zeroing \
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
VPUBLIC int setPartMG(MGparm *mgparm, Vpmg *pmg) {

    int j;
    double partMin[3], partMax[3];

    if (mgparm->type == 2) {
        for (j=0; j<3; j++) {
            partMin[j] = mgparm->center[j] + mgparm->partDisjCenterShift[j]
              - 0.5*mgparm->partDisjLength[j];
            partMax[j] = mgparm->center[j] + mgparm->partDisjCenterShift[j]
              + 0.5*mgparm->partDisjLength[j];
        }
        Vnm_print(0, "main:  Disj part lower corner = (%g, %g, %g)\n",
          partMin[0], partMin[1], partMin[2]);
        Vnm_print(0, "main:  Disj part upper corner = (%g, %g, %g)\n",
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
            Vnm_tprint( 1, "main:    Total electrostatic energy = \
%1.12E kJ/mol\n", Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*totEnergy));
        } else *totEnergy = 0;
    } else if (pbeparm->calcenergy == 2) {
        *nenergy = 1;
        *totEnergy = Vpmg_energy(pmg, extEnergy);
        *qfEnergy = Vpmg_qfEnergy(pmg, extEnergy);
        *qmEnergy = Vpmg_qmEnergy(pmg, extEnergy);
        *dielEnergy = Vpmg_dielEnergy(pmg, extEnergy);
        Vnm_tprint( 1, "main:    Total electrostatic energy = %1.12E \
kJ/mol\n", Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*totEnergy));
        Vnm_tprint( 1, "main:    Fixed charge energy = %g kJ/mol\n",
           0.5*Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*qfEnergy));
        Vnm_tprint( 1, "main:    Mobile charge energy = %g kJ/mol\n",
           Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*qmEnergy));
        Vnm_tprint( 1, "main:    Dielectric energy = %g kJ/mol\n",
           Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*dielEnergy));
        Vnm_tprint( 1, "main:    Per-atom energies:\n");
        alist = pmg->pbe->alist;
        for (i=0; i<Valist_getNumberAtoms(alist); i++) {
            atom = Valist_getAtom(alist, i); 
            tenergy = Vpmg_qfAtomEnergy(pmg, atom);
            Vnm_tprint( 1, "main:        Atom %d:  %1.12E kJ/mol\n", i,
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
// Args:     com         communications object
//           nosh        stores input file information
//           pbeparm     PBE parameters
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
                Vpmg_qfForce(pmg, qfForce, j);
                Vpmg_ibForce(pmg, ibForce, j);
                Vpmg_dbnpForce(pmg, dbForce, npForce, pbeparm->gamma, j);
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
        Vnm_tprint( 1, "main:    Net fixed charge force on molecule %d\n",
          pbeparm->molid);
        Vnm_tprint( 1, "           = (%4.3e, %4.3e, %4.3e) kJ/(mol A)\n",
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].qfForce[0],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].qfForce[1],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].qfForce[2]);
        Vnm_tprint( 1, "main:    Net ionic boundary force on molecule %d\n",
          pbeparm->molid);
        Vnm_tprint( 1, "           = (%4.3e, %4.3e, %4.3e) kJ/(mol A)\n",
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].ibForce[0],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].ibForce[1],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].ibForce[2]);
        Vnm_tprint( 1, "main:    Net dielectric boundary force on \
molecule %d\n", pbeparm->molid);
        Vnm_tprint( 1, "           = (%4.3e, %4.3e, %4.3e) kJ/(mol A)\n",
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].dbForce[0],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].dbForce[1],
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].dbForce[2]);
        Vnm_tprint( 1, "main:    Net apolar force on molecule %d\n",
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
                Vpmg_qfForce(pmg, (*atomForce)[j].qfForce, j);
                Vpmg_ibForce(pmg, (*atomForce)[j].ibForce, j);
                Vpmg_dbnpForce(pmg, (*atomForce)[j].dbForce,
                  (*atomForce)[j].npForce, pbeparm->gamma, j);
            } else {
                for (k=0; k<3; k++) {
                    (*atomForce)[j].qfForce[k] = 0;
                    (*atomForce)[j].ibForce[k] = 0;
                    (*atomForce)[j].dbForce[k] = 0;
                    (*atomForce)[j].npForce[k] = 0;
                }
            }
            Vnm_tprint( 1, "main:    Total force on atom %d, molecule %d \
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
            Vnm_tprint( 1, "main:    Fixed charge force on atom %d, \
molecule %d = (%4.3e, %4.3e, %4.3e) kJ/mol/A\n", j, pbeparm->molid,
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].qfForce[0],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].qfForce[1],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].qfForce[2]);
            Vnm_tprint( 1, "main:    Ionic boundary force on atom %d, \
molecule %d = (%4.3e, %4.3e, %4.3e) kJ/mol/A\n", j, pbeparm->molid,
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].ibForce[0],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].ibForce[1],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].ibForce[2]);
            Vnm_tprint( 1, "main:    Dielectric boundary force on atom \
%d, molecule %d = (%4.3e, %4.3e, %4.3e) kJ/mol/A\n", j, pbeparm->molid,
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].dbForce[0],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].dbForce[1],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].dbForce[2]);
            Vnm_tprint( 1, "main:    Apolar force on atom %d, molecule %d \
= (%4.3e, %4.3e, %4.3e) kJ/mol/A\n", j, pbeparm->molid,
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].npForce[0],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].npForce[1],
             Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[j].npForce[2]);
        }
    } else *nforce = 0;

    return 1;
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
    char outpath[72];
    char mxtype[3];

    if (nosh->bogus) return 1;

#ifdef HAVE_MPI_H
    sprintf(writematstem, "%s-PE%d", pbeparm->writematstem,
      rank);
#else
    sprintf(writematstem, "%s", pbeparm->writematstem);
#endif
    
    if (pbeparm->writemat == 1) {
        if (snprintf(outpath, 72, "%s.%s", writematstem, "mat") == -1) {
            Vnm_print(2, "main:    Matrix output path truncated to: %s\n", 
              outpath);
            Vnm_print(2, "main:    72-character limit exceeded!\n");
            Vnm_print(2, "main:    Skipping matrix I/O!\n");
            return 1;
        }
        mxtype[0] = 'R';
        mxtype[1] = 'S';
        mxtype[2] = 'A';
        /* Poisson operator only */
        if (pbeparm->writematflag == 0) {
            Vnm_tprint( 1, "main:    Writing Poisson operator matrix \
to %s...\n", outpath);

         /* Linearization of Poisson-Boltzmann operator around solution */
         } else if (pbeparm->writematflag == 1) {
            Vnm_tprint( 1, "main:    Writing linearization of full \
Poisson-Boltzmann operator matrix to %s...\n", outpath);

         } else {
             Vnm_tprint( 2, "main:    Bogus matrix specification\
(%d)!\n", pbeparm->writematflag);
             return 0;
         }

         Vnm_tprint(0, "main:    Printing operator...\n");
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
 
                Vnm_print(1, "main:    Writing charge distribution to ");
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
 
                Vnm_print(1, "main:    Writing potential to ");
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

                Vnm_print(1, "main:    Writing molecular accessibility to ");
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

                Vnm_print(1, "main:    Writing spline-based accessibility to ");
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

                Vnm_print(1, "main:    Writing van der Waals accessibility to ");
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

                Vnm_print(1, "main:    Writing ion accessibility to ");
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

                Vnm_print(1, "main:    Writing potential Laplacian to ");
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

                Vnm_print(1, "main:    Writing energy density to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                Vpmg_fillArray(pmg, pmg->rwork, VDT_EDENS, 0.0);
                sprintf(title, 
                  "ENERGY DENSITY (kT/e/A)^2");
                break;

            case VDT_NDENS:

                Vnm_print(1, "main:    Writing number density to ");
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

                Vnm_print(1, "main:    Writing charge density to ");
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

                Vnm_print(1, "main:    Writing x-shifted dielectric map to ");
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

                Vnm_print(1, "main:    Writing y-shifted dielectric map to ");
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

                Vnm_print(1, "main:    Writing z-shifted dielectric map to ");
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

                Vnm_print(1, "main:    Writing kappa map to ");
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

                Vnm_print(2, "Invalid data type for writing!\n");
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
                Vnm_print(1, "%s\n", outpath);
                grid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin,
                  pmg->rwork);
                Vgrid_writeDX(grid, "FILE", "ASC", VNULL, outpath, title,
                  pmg->pvec);
                Vgrid_dtor(&grid);
                break;

            case VDF_AVS:
                sprintf(outpath, "%s.%s", writestem, "ucd");
                Vnm_print(1, "%s\n", outpath);
                Vnm_print(2, "main:  Sorry, AVS format isn't supported for \
uniform meshes yet!\n");
                break;

            case VDF_UHBD:
                sprintf(outpath, "%s.%s", writestem, "grd");
                Vnm_print(1, "%s\n", outpath);
                grid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin,
                  pmg->rwork);
                Vgrid_writeUHBD(grid, "FILE", "ASC", VNULL, outpath, title,
                  pmg->pvec);
                Vgrid_dtor(&grid);
                break;

            default:
                Vnm_print(2, "main:  Bogus data format (%d)!\n", 
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
    double ltenergy, gtenergy;

    Vnm_tprint( 1, "main:  print energy %d ", nosh->printcalc[i][0]);
    for (j=1; j<nosh->printnarg[i]; j++) {
        if (nosh->printop[i][j-1] == 0)
          Vnm_print(1, "+ ");
        else if (nosh->printop[i][j-1] == 1)
          Vnm_print(1, "- ");
        else {
            Vnm_tprint( 2, "main:  Undefined PRINT operation!\n");
            return 0;
        }
        Vnm_print(1, "%d ", nosh->printcalc[i][j]);
    }
    Vnm_print(1, "end\n");
    calcid = nosh->elec2calc[nosh->printcalc[i][0]-1];
    if (nosh->calc[calcid].pbeparm->calcenergy != 0) {
        ltenergy = Vunit_kb * (1e-3) * Vunit_Na *
          nosh->calc[calcid].pbeparm->temp * totEnergy[calcid];
    } else {
        Vnm_tprint( 2, "main:    Didn't calculate energy in Calculation \
#%d\n", calcid+1);
        return 0;
    }
    for (j=1; j<nosh->printnarg[i]; j++) {
        calcid = nosh->elec2calc[nosh->printcalc[i][j]-1];
        if (nosh->calc[calcid].pbeparm->calcenergy != 0) {
            if (nosh->printop[i][j-1] == 0)
              ltenergy = ltenergy + Vunit_kb * (1e-3) * Vunit_Na *
                nosh->calc[calcid].pbeparm->temp * totEnergy[calcid];
            else if (nosh->printop[i][j-1] == 1)
              ltenergy = ltenergy - Vunit_kb * (1e-3) * Vunit_Na *
                nosh->calc[calcid].pbeparm->temp * totEnergy[calcid];
        } else {
            Vnm_tprint( 2, "main:    Didn't calculate energy in \
Calculation #%d\n", calcid+1);
            return 0;
        }
    }

    Vnm_tprint( 1, "main:    Local answer (PE %d) = %1.12E kJ/mol\n", 
      Vcom_rank(com), ltenergy);
    Vnm_tprint( 0, "printEnergy:  Performing global reduction (sum)\n");
    Vcom_reduce(com, &ltenergy, &gtenergy, 1, 2, 0);
    Vnm_tprint( 1, "main:    Global answer = %1.12E kJ/mol\n", gtenergy);

    return 1;

}

