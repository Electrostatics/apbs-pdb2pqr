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

#include "apbscfg.h"
#include "maloc/maloc.h"  
#ifdef HAVE_MC_H
#  include "mc/mc.h"  
#endif
#ifdef HAVE_MCX_H
#  include "mcx/mcx.h"  
#endif

#include "apbs/apbs.h"  
#include "apbs/vhal.h"  
#include "apbs/nosh.h"  
#include "apbs/vgrid.h"  
#include "apbs/mgparm.h"  
#include "apbs/pbeparm.h"  
#include "apbs/femparm.h"  

#include "routines.h"

VEMBED(rcsid="$Id$")

VPUBLIC void startVio() { Vio_start(); }

VPUBLIC int loadMolecules(NOsh *nosh, Valist *alist[NOSH_MAXMOL]) {
    
    int i, j, rc;
    double q; 
    Vatom *atom = VNULL;
    Vparam *param = VNULL;

    Vnm_tprint( 1, "Got PQR paths for %d molecules\n", nosh->nmol);
    if (nosh->nmol <= 0) {
       Vnm_tprint(2, "You didn't specify any molecules (correctly)!\n");
       Vnm_tprint(2, "Bailing out!\n");
       return 0;
    }

    if (nosh->gotparm) {
        param = Vparam_ctor();
        switch (nosh->parmfmt) {
            case NPF_FLAT:
                 Vnm_tprint( 1, "Reading parameter data from %s.\n",
                   nosh->parmpath);
                if (Vparam_readFlatFile(param, "FILE", "ASC", VNULL, 
                  nosh->parmpath) != 1) {
                    Vnm_tprint(2, "NOsh:  Error reading parameter\
 file (%s)!\n", nosh->parmpath);
                    return 0;
                }
                break;
            default:
                Vnm_tprint(2, "NOsh:  Error! Undefined parameter file \
type (%d)!\n", nosh->parmfmt);
                return 0;
        } /* switch parmfmt */
    }

    for (i=0; i<nosh->nmol; i++) {
        alist[i] = Valist_ctor();
        switch (nosh->molfmt[i]) {
            case NMF_PQR:
                Vnm_tprint( 1, "Reading PQR-format atom data from %s.\n",
                  nosh->molpath[i]);
                rc = Valist_readPQR(alist[i], "FILE", "ASC", VNULL,
                  nosh->molpath[i]);
                break;
            case NMF_PDB:
                /* Load parameters */
                if (!nosh->gotparm) {
                    Vnm_tprint(2, "NOsh:  Error!  Can't read PDB without \
specifying PARM file!\n");
                    return 0;
                }
                Vnm_tprint( 1, "Reading PDB-format atom data from %s.\n",
                  nosh->molpath[i]);
                rc = Valist_readPDB(alist[i], param, "FILE", "ASC", VNULL,
                  nosh->molpath[i]);
                break;
            default:
                Vnm_tprint(2, "NOsh:  Error!  Undefined molecule file type \
(%d)!\n", nosh->molfmt[i]);
                return 0;
        } /* switch molfmt */

        if (rc != 1) {
            Vnm_tprint( 2, "Error while reading molecule from %s\n",
              nosh->molpath[i]);
            return 0;
        }

        Vnm_tprint( 1, "  %d atoms\n", Valist_getNumberAtoms(alist[i]));
        Vnm_tprint( 1, "  Centered at (%4.3e, %4.3e, %4.3e)\n",
          alist[i]->center[0], alist[i]->center[1], alist[i]->center[2]);
        Vnm_tprint( 1, "  Net charge %3.2e e\n", alist[i]->charge);        

        /* Check for uncharged molecule */
        q = 0;
        for (j=0; j<Valist_getNumberAtoms(alist[i]); j++) {
            atom = Valist_getAtom(alist[i], j);
            q += VSQR(Vatom_getCharge(atom));
        }
        if (q < (1e-6)) {
            Vnm_tprint(2, "Molecule #%d is uncharged!\n");
            Vnm_tprint(2, "Sum square charge = %g\n", q);
            Vnm_tprint(2, "Bailing out!\n");
        }
    }

    if (nosh->gotparm) Vparam_dtor(&param);

    return 1;

}

VPUBLIC void killMolecules(NOsh *nosh, Valist *alist[NOSH_MAXMOL]) {
    
    int i;

#ifndef VAPBSQUIET
    Vnm_tprint( 1, "Destroying %d molecules\n", nosh->nmol);
#endif

    for (i=0; i<nosh->nmol; i++) Valist_dtor(&(alist[i]));

}

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
        switch (nosh->dielfmt[i]) {
            case VDF_DX:
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
                Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);
                break;
            case VDF_UHBD:
                Vnm_tprint( 2, "UHBD input not supported yet!\n");
                return 0;
            case VDF_AVS:
                Vnm_tprint( 2, "AVS input not supported yet!\n");
                return 0;
            default:
                Vnm_tprint( 2, "Invalid data format (%d)!\n", 
                  nosh->dielfmt[i]);
                return 0;
        }
        Vnm_tprint( 1, "Reading y-shifted dielectric map data from \
%s:\n", nosh->dielYpath[i]);
        dielYMap[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        switch (nosh->dielfmt[i]) {
            case VDF_DX:
                if (Vgrid_readDX(dielYMap[i], "FILE", "ASC", VNULL, 
                  nosh->dielYpath[i]) != 1) {
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
                Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);
                break;
            case VDF_UHBD:
                Vnm_tprint( 2, "UHBD input not supported yet!\n");
                return 0;
            case VDF_AVS:
                Vnm_tprint( 2, "AVS input not supported yet!\n");
                return 0;
            default:
                Vnm_tprint( 2, "Invalid data format (%d)!\n", 
                  nosh->dielfmt[i]);
                return 0;
        }
        Vnm_tprint( 1, "Reading z-shifted dielectric map data from \
%s:\n", nosh->dielZpath[i]);
        dielZMap[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        switch (nosh->dielfmt[i]) {
            case VDF_DX:
                if (Vgrid_readDX(dielZMap[i], "FILE", "ASC", VNULL, 
                  nosh->dielZpath[i]) != 1) {
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
                for (ii=0; ii<(nx*ny*nz); ii++) sum += (dielZMap[i]->data[ii]);
                sum = sum*hx*hy*hzed;
                Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);
                break;
            case VDF_UHBD:
                Vnm_tprint( 2, "UHBD input not supported yet!\n");
                return 0;
            case VDF_AVS:
                Vnm_tprint( 2, "AVS input not supported yet!\n");
                return 0;
            default:
                Vnm_tprint( 2, "Invalid data format (%d)!\n", 
                  nosh->dielfmt[i]);
                return 0;
        }
    }

    return 1;

}

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
        switch (nosh->kappafmt[i]) {
            case VDF_DX:
                if (Vgrid_readDX(map[i], "FILE", "ASC", VNULL, 
                  nosh->kappapath[i]) != 1) {
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
                break;
            case VDF_UHBD:
                Vnm_tprint( 2, "UHBD input not supported yet!\n");
                return 0;
            case VDF_AVS:
                Vnm_tprint( 2, "AVS input not supported yet!\n");
                return 0;
            default:
                Vnm_tprint( 2, "Invalid data format (%d)!\n", 
                  nosh->kappafmt[i]);
                return 0;
        }
    }

    return 1;

}

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
        switch (nosh->chargefmt[i]) {
            case VDF_DX:
                if (Vgrid_readDX(map[i], "FILE", "ASC", VNULL, 
                  nosh->chargepath[i]) != 1) {
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
                break;
            case VDF_UHBD:
                Vnm_tprint( 2, "UHBD input not supported yet!\n");
                return 0;
            case VDF_AVS:
                Vnm_tprint( 2, "AVS input not supported yet!\n");
                return 0;
            default:
                Vnm_tprint( 2, "Invalid data format (%d)!\n", 
                  nosh->kappafmt[i]);
                return 0;
        }
    }

    return 1;

}

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

VPUBLIC void printPBEPARM(PBEparm *pbeparm) {
    
    int i;
    double ionstr = 0.0;

    for (i=0; i<pbeparm->nion; i++)
      ionstr += 0.5*(VSQR(pbeparm->ionq[i])*pbeparm->ionc[i]);

    Vnm_tprint( 1, "  Molecule ID: %d\n", pbeparm->molid);
    switch (pbeparm->pbetype) {
        case PBE_NPBE:
            Vnm_tprint( 1, "  Nonlinear traditional PBE\n");
            break;
        case PBE_LPBE:
            Vnm_tprint( 1, "  Linearized traditional PBE\n");
            break;
        case PBE_NRPBE:
            Vnm_tprint( 1, "  Nonlinear regularized PBE\n");
            break;
        case PBE_LRPBE:
            Vnm_tprint( 1, "  Linearized regularized PBE\n");
            break;
        default:
            Vnm_tprint(2, "  Unknown PBE type (%d)!\n", pbeparm->pbetype);
            break;
    }
    if (pbeparm->bcfl == BCFL_ZERO) {
        Vnm_tprint( 1, "  Zero boundary conditions\n");
    } else if (pbeparm->bcfl == BCFL_SDH) {
        Vnm_tprint( 1, "  Single Debye-Huckel sphere boundary \
conditions\n");
    } else if (pbeparm->bcfl == BCFL_MDH) {
        Vnm_tprint( 1, "  Multiple Debye-Huckel sphere boundary \
conditions\n");
    } else if (pbeparm->bcfl == BCFL_FOCUS) {
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

VPUBLIC void printMGPARM(MGparm *mgparm, double realCenter[3]) {

    switch (mgparm->chgm) {
        case 0:
            Vnm_tprint(1, "  Using linear spline charge discretization.\n");
            break;
        case 1:
            Vnm_tprint(1, "  Using cubic spline charge discretization.\n");
            break;
        default:
            break;
    }
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

VPUBLIC int initMG(int i, NOsh *nosh, MGparm *mgparm, 
  PBEparm *pbeparm, double realCenter[3], Vpbe *pbe[NOSH_MAXCALC], 
  Valist *alist[NOSH_MAXMOL], Vgrid *dielXMap[NOSH_MAXMOL], 
  Vgrid *dielYMap[NOSH_MAXMOL], Vgrid *dielZMap[NOSH_MAXMOL],
  Vgrid *kappaMap[NOSH_MAXMOL], Vgrid *chargeMap[NOSH_MAXMOL],
  Vpmgp *pmgp[NOSH_MAXCALC], Vpmg *pmg[NOSH_MAXCALC]) {
    
    int j, bytesTotal, highWater, imol, focusFlag;
    double sparm, iparm;
    Vgrid *theDielXMap, *theDielYMap, *theDielZMap, *theKappaMap, *theChargeMap;

    Vnm_tstart(27, "Setup timer");

    /* Fix mesh center for "GCENT MOL #" types of declarations. */
    if (mgparm->cmeth == MCM_MOL) {
        Vnm_tprint(0, "Fixing grid center based on molecule...\n");
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
    Vnm_tprint(0, "Fixing grid center...\n");
    if (mgparm->type == 2) {
        for (j=0; j<3; j++) realCenter[j] = mgparm->center[j]
          + mgparm->partOlapCenterShift[j];
    } else {
        for (j=0; j<3; j++) realCenter[j] = mgparm->center[j];
    }

    /* Set up PBE object */
    Vnm_tprint(0, "Setting up PBE object...\n");
    if (pbeparm->srfm == VSM_SPLINE) sparm = pbeparm->swin;
    else sparm = pbeparm->srad;
    if (pbeparm->nion > 0) iparm = pbeparm->ionr[0];
    else iparm = 0.0;
	if (pbeparm->bcfl == BCFL_FOCUS) {
	  if (i == 0) {
            Vnm_tprint( 2, "Can't focus first calculation!\n");
            return 0;
        }
	  focusFlag = 1;
	} else focusFlag = 0;
	
	pbe[i] = Vpbe_ctor(alist[pbeparm->molid-1], pbeparm->nion,
		  pbeparm->ionc, pbeparm->ionr, pbeparm->ionq, pbeparm->temp,
		  pbeparm->gamma, pbeparm->pdie, pbeparm->sdie, sparm, focusFlag);

    /* Set up PDE object */
    Vnm_tprint(0, "Setting up PDE object...\n");
    switch (pbeparm->pbetype) {
        case PBE_NPBE:
            pmgp[i] = Vpmgp_ctor(mgparm->dime[0], mgparm->dime[1],
              mgparm->dime[2], mgparm->nlev, mgparm->grid[0], mgparm->grid[1],
              mgparm->grid[2], 1);
            break;
        case PBE_LPBE:
            pmgp[i] = Vpmgp_ctor(mgparm->dime[0], mgparm->dime[1],
              mgparm->dime[2], mgparm->nlev, mgparm->grid[0], mgparm->grid[1],
              mgparm->grid[2], 0);
            break;
        case PBE_LRPBE:
            Vnm_tprint(2, "Sorry, LRPBE isn't supported with the MG solver!\n");
            return 0;
            break;
        case PBE_NRPBE:
            Vnm_tprint(2, "Sorry, NRPBE isn't supported with the MG solver!\n");
            return 0;
            break;
        default:
            Vnm_tprint(2, "Error!  Unknown PBE type (%d)!\n", pbeparm->pbetype);
            return 0;
    }
    Vnm_tprint(0, "Setting PDE center to local center...\n");
    pmgp[i]->bcfl = pbeparm->bcfl;
    pmgp[i]->xcent = realCenter[0];
    pmgp[i]->ycent = realCenter[1];
    pmgp[i]->zcent = realCenter[2];
    if (pbeparm->bcfl == BCFL_FOCUS) {
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
    if (pbeparm->useDielMap) {
        if ((pbeparm->dielMapID-1) < nosh->ndiel) {
            theDielXMap = dielXMap[pbeparm->dielMapID-1];
        } else {
            Vnm_print(2, "Error!  %d is not a valid dielectric map ID!\n", 
                    pbeparm->dielMapID);
            return 0;
        }
    } else theDielXMap = VNULL;
    if (pbeparm->useDielMap) {
        if ((pbeparm->dielMapID-1) < nosh->ndiel) {
            theDielYMap = dielYMap[pbeparm->dielMapID-1];
        } else {
            Vnm_print(2, "Error!  %d is not a valid dielectric map ID!\n",
                    pbeparm->dielMapID);
            return 0;
        }
    } else theDielYMap = VNULL;
    if (pbeparm->useDielMap) {
        if ((pbeparm->dielMapID-1) < nosh->ndiel) {
            theDielZMap = dielZMap[pbeparm->dielMapID-1];
        } else {
            Vnm_print(2, "Error!  %d is not a valid dielectric map ID!\n",
                    pbeparm->dielMapID);
            return 0;
        }
    } else theDielZMap = VNULL;
    if (pbeparm->useKappaMap) {
        if ((pbeparm->kappaMapID-1) < nosh->nkappa) {
            theKappaMap = kappaMap[pbeparm->kappaMapID-1];
        } else {
            Vnm_print(2, "Error!  %d is not a valid kappa map ID!\n",
                    pbeparm->kappaMapID);
            return 0;
        }
    } else theKappaMap = VNULL;
    if (pbeparm->useChargeMap) {
        if ((pbeparm->chargeMapID-1) < nosh->ncharge) {
            theChargeMap = chargeMap[pbeparm->chargeMapID-1];
        } else {
            Vnm_print(2, "Error!  %d is not a valid charge map ID!\n",
                    pbeparm->chargeMapID);
            return 0;
        }
    } else theChargeMap = VNULL;
    Vpmg_fillco(pmg[i], 
      pbeparm->srfm, pbeparm->swin, mgparm->chgm,
      pbeparm->useDielMap, theDielXMap,
      pbeparm->useDielMap, theDielYMap,
      pbeparm->useDielMap, theDielZMap,
      pbeparm->useKappaMap, theKappaMap,
      pbeparm->useChargeMap, theChargeMap);

    /* Print a few derived parameters */
#ifndef VAPBSQUIET
    Vnm_tprint(1, "  Debye length:  %g A\n", Vpbe_getDeblen(pbe[i]));
#endif

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

VPUBLIC void killMG(NOsh *nosh, Vpbe *pbe[NOSH_MAXCALC], 
  Vpmgp *pmgp[NOSH_MAXCALC], Vpmg *pmg[NOSH_MAXCALC]) {
    
#ifndef VAPBSQUIET
    Vnm_tprint(1, "Destroying multigrid structures.\n");
#endif

    Vpbe_dtor(&(pbe[nosh->ncalc-1]));
    Vpmg_dtor(&(pmg[nosh->ncalc-1]));
    Vpmgp_dtor(&(pmgp[nosh->ncalc-1]));

}

VPUBLIC int solveMG(NOsh *nosh, Vpmg *pmg, MGparm_CalcType type) {

    int nx, ny, nz, i;

   
    if (nosh != VNULL) {
        if (nosh->bogus) return 1;
    }

    Vnm_tstart(28, "Solver timer");


    if (type != MCT_DUM) {
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

VPUBLIC int forceMG(Vmem *mem, NOsh *nosh, PBEparm *pbeparm, MGparm *mgparm,
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
                Vpmg_qfForce(pmg, qfForce, j, mgparm->chgm);
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
                Vpmg_qfForce(pmg, (*atomForce)[j].qfForce, j, mgparm->chgm);
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

VPUBLIC void killEnergy() { 

#ifndef VAPBSQUIET
    Vnm_tprint(1, "No energy arrays to destroy.\n"); 
#endif

}

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

VPUBLIC int printEnergy(Vcom *com, NOsh *nosh, double totEnergy[NOSH_MAXCALC], 
  int i) {

    int j, calcid;
    double ltenergy, gtenergy, scalar;
    
    if (Vstring_strcasecmp(nosh->elecname[nosh->printcalc[i][0]], "") == 0){
      Vnm_tprint( 1, "print energy %d ", nosh->printcalc[i][0]);
    } else {
      Vnm_tprint( 1, "print energy %d (%s) ", nosh->printcalc[i][0], nosh->elecname[nosh->printcalc[i][0]]);
    }
    for (j=1; j<nosh->printnarg[i]; j++) {
        if (nosh->printop[i][j-1] == 0)
          Vnm_tprint(1, "+ ");
        else if (nosh->printop[i][j-1] == 1)
          Vnm_tprint(1, "- ");
        else {
            Vnm_tprint( 2, "Undefined PRINT operation!\n");
            return 0;
        }
        if (Vstring_strcasecmp(nosh->elecname[j+1], "") == 0){
          Vnm_tprint( 1, "%d ", nosh->printcalc[i][j]);
        } else {
          Vnm_tprint( 1, "%d (%s) ", nosh->printcalc[i][j], nosh->elecname[nosh->printcalc[i][j]]);
        }
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

VPUBLIC int printForce(Vcom *com, NOsh *nosh, int nforce[NOSH_MAXCALC], 
  AtomForce *atomForce[NOSH_MAXCALC], int i) {

    int ipr, ifr, ivc, calcid, refnforce, refcalcforce;
    double temp, scalar;
    AtomForce *lforce, *gforce, *aforce;

    if (Vstring_strcasecmp(nosh->elecname[nosh->printcalc[i][0]], "") == 0){
      Vnm_tprint( 1, "print force %d ", nosh->printcalc[i][0]);
    } else {
      Vnm_tprint( 1, "print force %d (%s) ", nosh->printcalc[i][0], nosh->elecname[nosh->printcalc[i][0]]);
    }
    for (ipr=1; ipr<nosh->printnarg[i]; ipr++) {
        if (nosh->printop[i][ipr-1] == 0)
          Vnm_tprint(1, "+ ");
        else if (nosh->printop[i][ipr-1] == 1)
          Vnm_tprint(1, "- ");
        else {
            Vnm_tprint( 2, "Undefined PRINT operation!\n");
            return 0;
        }
        if (Vstring_strcasecmp(nosh->elecname[ipr+1], "") == 0){
          Vnm_tprint( 1, "%d ", nosh->printcalc[i][ipr]);
        } else {
          Vnm_tprint( 1, "%d (%s) ", nosh->printcalc[i][ipr], nosh->elecname[nosh->printcalc[i][ipr]]);
        }
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

#ifdef HAVE_MC_H

VPUBLIC int initFE(int icalc, NOsh *nosh, FEMparm *feparm, PBEparm *pbeparm, 
  Vpbe *pbe[NOSH_MAXCALC], Valist *alist[NOSH_MAXMOL], 
  Vfetk *fetk[NOSH_MAXCALC]) {
    
    int j, bytesTotal, highWater, theMol, focusFlag;
    double sparm, iparm, center[3];

    Vnm_tstart(27, "Setup timer");

    /* Print some warning messages */
    if (pbeparm->useDielMap)  Vnm_tprint(2, "FEM ignoring dielectric map!\n");
    if (pbeparm->useKappaMap)  Vnm_tprint(2, "FEM ignoring kappa map!\n");
    if (pbeparm->useChargeMap)  Vnm_tprint(2, "FEM ignoring charge map!\n");

    /* Fix mesh center for "GCENT MOL #" types of declarations. */
    Vnm_tprint(0, "Re-centering mesh...\n");
    theMol = pbeparm->molid-1;
    for (j=0; j<3; j++) {
        if (theMol < nosh->nmol) {
            center[j] = (alist[theMol])->center[j];
        } else{ 
            Vnm_tprint(2, "ERROR!  Bogus molecule number (%d)!\n", 
              (theMol+1));
            return 0;
        }
    }
    
    /* Set up PBE object */
    Vnm_tprint(0, "Setting up PBE object...\n");
    if (pbeparm->srfm == VSM_SPLINE) sparm = pbeparm->swin;
    else sparm = pbeparm->srad;
    if (pbeparm->nion > 0) iparm = pbeparm->ionr[0];
    else iparm = 0.0;
	focusFlag = 0;
    pbe[icalc] = Vpbe_ctor(alist[theMol], pbeparm->nion,
      pbeparm->ionc, pbeparm->ionr, pbeparm->ionq, pbeparm->temp,
      pbeparm->gamma, pbeparm->pdie, pbeparm->sdie, sparm, focusFlag);

    /* Print a few derived parameters */
    Vnm_tprint(1, "  Debye length:  %g A\n", Vpbe_getDeblen(pbe[icalc]));

    /* Set up FEtk objects */
    Vnm_tprint(0, "Setting up FEtk object...\n");
    fetk[icalc] = Vfetk_ctor(pbe[icalc], PBE_NRPBE);
    Vfetk_setParameters(fetk[icalc], pbeparm, feparm);

    /* Build mesh */
    Vnm_tprint(0, "Setting up mesh...\n");
    Vfetk_genCube(fetk[icalc], alist[theMol]->center, feparm->domainLength);
    /* Uniformly refine the mesh a bit */
    for (j=0; j<2; j++) {
        AM_markRefine(fetk[icalc]->am, 0, -1, 0, 0);
        AM_refine(fetk[icalc]->am, 2, USEHB);
        Vnm_redirect(0);
        Gem_shapeChk(fetk[icalc]->gm);
        Vnm_redirect(1);
    }

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

VPUBLIC void printFEPARM(int icalc, NOsh *nosh, FEMparm *feparm, 
  Vfetk *fetk[NOSH_MAXCALC]) {

    Vnm_tprint(1, "  Domain size:  %g A x %g A x %g A\n", 
      feparm->domainLength[0], feparm->domainLength[1],
      feparm->domainLength[2]);
    switch(feparm->ekey) {
        case FET_SIMP:
            Vnm_tprint(1, "  Per-simplex error tolerance:  %g\n", feparm->etol);
            break;
        case FET_GLOB:
            Vnm_tprint(1, "  Global error tolerance:  %g\n", feparm->etol);
            break;
        case FET_FRAC:
            Vnm_tprint(1, "  Fraction of simps to refine:  %g\n", feparm->etol);
            break;
        default:
            Vnm_tprint(2, "Invalid ekey (%d)!\n", feparm->ekey);
            VASSERT(0);
            break;
    }
    switch(feparm->akeyPRE) {
        case FRT_UNIF:
            Vnm_tprint(1, "  Uniform pre-solve refinement.\n");
            break;
        case FRT_GEOM:
            Vnm_tprint(1, "  Geometry-based pre-solve refinement.\n");
            break;
        case FRT_RESI:
            Vnm_tprint(1, "  Residual-based pre-solve refinement.\n");
            Vnm_tprint(2, "What?  You can't do a posteriori error estimation \
before you solve!\n");
            VASSERT(0);
            break;
        case FRT_DUAL:
            Vnm_tprint(1, "  Dual-based pre-solve refinement.\n");
            Vnm_tprint(2, "What?  You can't do a posteriori error estimation \
before you solve!\n");
            VASSERT(0);
            break;
        case FRT_LOCA:
            Vnm_tprint(1, "  Local-based pre-solve refinement.\n");
            Vnm_tprint(2, "What?  You can't do a posteriori error estimation \
before you solve!\n");
            VASSERT(0);
            break;
        default:
            Vnm_tprint(2, "Invalid akeyPRE (%d)!\n", feparm->akeyPRE);
            VASSERT(0);
            break;
    }
    switch(feparm->akeySOLVE) {
        case FRT_UNIF:
            Vnm_tprint(1, "  Uniform a posteriori refinement.\n");
            break;
        case FRT_GEOM:
            Vnm_tprint(1, "  Geometry-based a posteriori refinement.\n");
            break;
        case FRT_RESI:
            Vnm_tprint(1, "  Residual-based a posteriori refinement.\n");
            break;
        case FRT_DUAL:
            Vnm_tprint(1, "  Dual-based a posteriori refinement.\n");
            break;
        case FRT_LOCA:
            Vnm_tprint(1, "  Local-based a posteriori refinement.\n");
            break;
        default:
            Vnm_tprint(2, "Invalid akeySOLVE (%d)!\n", feparm->akeySOLVE);
            break;
    }
    Vnm_tprint(1, "  Refinement of initial mesh to ~%d vertices\n", 
      feparm->targetNum);
    Vnm_tprint(1, "  Geometry-based refinment lower bound:  %g A\n",
      feparm->targetRes);
    Vnm_tprint(1, "  Maximum number of solve-estimate-refine cycles:  %d\n",
      feparm->maxsolve);
    Vnm_tprint(1, "  Maximum number of vertices in mesh:  %d\n",
      feparm->maxvert);

    /* FOLLOWING IS SOLVER-RELATED; BAIL IF NOT SOLVING */
    if (nosh->bogus)  return;
    if (USEHB) {
        Vnm_tprint(1, "  HB linear solver:  AM_hPcg\n");
    } else {
        Vnm_tprint(1, "  Non-HB linear solver:  ");
        switch (fetk[icalc]->lkey) {
            case VLT_SLU:
                Vnm_print(1, "SLU direct\n");
                break;
            case VLT_MG:
                Vnm_print(1, "multigrid\n");
                break;
            case VLT_CG:
                Vnm_print(1, "conjugate gradient\n");
                break;
            case VLT_BCG:
                Vnm_print(1, "BiCGStab\n");
                break;
            default:
                Vnm_print(1, "???\n");
                break;
        }
    }
    Vnm_tprint(1, "  Linear solver tol.:  %g\n", fetk[icalc]->ltol);
    Vnm_tprint(1, "  Linear solver max. iters.:  %d\n", fetk[icalc]->lmax);
    Vnm_tprint(1, "  Linear solver preconditioner:  ");
    switch (fetk[icalc]->lprec) {
        case VPT_IDEN:
            Vnm_print(1, "identity\n");
            break;
        case VPT_DIAG:
            Vnm_print(1, "diagonal\n");
            break;
        case VPT_MG:
            Vnm_print(1, "multigrid\n");
            break;
        default:
            Vnm_print(1, "???\n");
            break;
    }
    Vnm_tprint(1, "  Nonlinear solver:  ");
    switch (fetk[icalc]->nkey) {
        case VNT_NEW:
            Vnm_print(1, "newton\n");
            break;
        case VNT_INC:
            Vnm_print(1, "incremental\n");
            break;
        case VNT_ARC:
            Vnm_print(1, "pseudo-arclength\n");
            break;
        default:
            Vnm_print(1, "??? ");
            break;
    }
    Vnm_tprint(1, "  Nonlinear solver tol.:  %g\n", fetk[icalc]->ntol);
    Vnm_tprint(1, "  Nonlinear solver max. iters.:  %d\n", fetk[icalc]->nmax);
    Vnm_tprint(1, "     Initial guess:  ");
    switch (fetk[icalc]->gues) {
        case VGT_ZERO:
            Vnm_tprint(1, "zero\n");
            break;
        case VGT_DIRI:
            Vnm_tprint(1, "boundary function\n");
            break;
        case VGT_PREV:
            Vnm_tprint(1, "interpolated previous solution\n");
            break;
        default:
            Vnm_tprint(1, "???\n");
            break;
    }

}

VPUBLIC int partFE(int icalc, NOsh *nosh, FEMparm *feparm, 
  Vfetk *fetk[NOSH_MAXCALC]) {

    Vfetk_setAtomColors(fetk[icalc]);
    return 1;
}

VPUBLIC int preRefineFE(int icalc, NOsh *nosh, FEMparm *feparm, 
  Vfetk *fetk[NOSH_MAXCALC]) {

    int nverts, marked;

    switch(feparm->akeyPRE) {
        case FRT_UNIF:
            Vnm_tprint(1, "  Commencing uniform refinement to %d verts.\n",
              feparm->targetNum);
            break;
        case FRT_GEOM:
            Vnm_tprint(1, "  Commencing geometry-based refinement to %d \
verts or %g A resolution.\n", feparm->targetNum, feparm->targetRes);
            break;
        case FRT_RESI:
            VASSERT(0);
            break;
        case FRT_DUAL:
            Vnm_tprint(2, "What?  You can't do a posteriori error estimation \
before you solve!\n");
            VASSERT(0);
            break;
        case FRT_LOCA:
            VASSERT(0);
            break;
        default:
            VASSERT(0);
            break;
    }

    Vnm_tprint(1, "  Initial mesh has %d vertices\n", 
      Gem_numVV(fetk[icalc]->gm));
    while (1) {
        nverts = Gem_numVV(fetk[icalc]->gm);
        if (nverts > feparm->targetNum) {
            Vnm_tprint(1, "  Hit vertex number limit.\n");
            break;
        }
        Vnm_print(1, "DEBUG -- akeyPRE = %d\n", feparm->akeyPRE);
        marked = AM_markRefine(fetk[icalc]->am, feparm->akeyPRE, -1, 
          feparm->ekey, feparm->etol);
        if (marked == 0) {
            Vnm_tprint(1, "  Marked 0 simps; hit error/size tolerance.\n");
            break;
        }
        Vnm_tprint(1, "    Have %d verts, marked %d.  Refining...\n", nverts,
          marked);
        AM_refine(fetk[icalc]->am, 0, USEHB);
    }
    nverts = Gem_numVV(fetk[icalc]->gm);
    Vnm_tprint(1, "  Done refining; have %d verts.\n", nverts);

    return 1;
}

VPUBLIC int solveFE(int icalc, NOsh *nosh, PBEparm *pbeparm, FEMparm *feparm, 
  Vfetk *fetk[NOSH_MAXCALC]) {

    int lkeyHB = 3;  /**<  AM_hPcg */

    if ((pbeparm->pbetype==PBE_NPBE)||(pbeparm->pbetype == PBE_NRPBE)) {
        AM_nSolve(fetk[icalc]->am, fetk[icalc]->nkey, fetk[icalc]->nmax, 
          fetk[icalc]->ntol, fetk[icalc]->lkey, fetk[icalc]->lmax, 
          fetk[icalc]->ltol, fetk[icalc]->lprec, fetk[icalc]->gues, 
          fetk[icalc]->pjac);
    } else if ((pbeparm->pbetype==PBE_LPBE)||(pbeparm->pbetype==PBE_LRPBE)) {
        if (USEHB) {
            AM_hlSolve(fetk[icalc]->am, 0, lkeyHB, fetk[icalc]->lmax, 
              fetk[icalc]->ltol, fetk[icalc]->gues, fetk[icalc]->pjac);
        } else {
            AM_lSolve(fetk[icalc]->am, 0, fetk[icalc]->lkey, fetk[icalc]->lmax, 
              fetk[icalc]->ltol, fetk[icalc]->lprec, fetk[icalc]->gues, 
              fetk[icalc]->pjac);
        }
    }

    return 1;
}

VPUBLIC int energyFE(NOsh *nosh, int icalc, Vfetk *fetk[NOSH_MAXCALC], 
  int *nenergy, double *totEnergy, double *qfEnergy, double *qmEnergy,
  double *dielEnergy) {

    double tenergy;
    FEMparm *feparm;
    PBEparm *pbeparm;

    feparm = nosh->calc[icalc].femparm;
    pbeparm = nosh->calc[icalc].pbeparm;

    *nenergy = 1;

    /* Some processors don't count */
    if (nosh->bogus == 0) {
        if ((pbeparm->pbetype==PBE_NPBE)||(pbeparm->pbetype==PBE_NRPBE)) {
            *totEnergy = Vfetk_energy(fetk[icalc], -1, 1);
        } else if ((pbeparm->pbetype==PBE_LPBE)||(pbeparm->pbetype==PBE_LRPBE)) {
            *totEnergy = Vfetk_energy(fetk[icalc], -1, 0);
        } else VASSERT(0);

#ifndef VAPBSQUIET
        Vnm_tprint(1, "      Total electrostatic energy = %1.12E kJ/mol\n", 
          Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*totEnergy));
        fflush(stdout);
#endif
    } else *totEnergy = 0;

    if (pbeparm->calcenergy == 2) {

        Vnm_tprint(2, "Error!  Verbose energy evaluation not available for FEM yet!\n");
        Vnm_tprint(2, "E-mail baker@biochem.wustl.edu if you want this.\n");
        *qfEnergy = 0;
        *qmEnergy = 0;
        *dielEnergy = 0;

    } else *nenergy = 0;

    return 1;
}

VPUBLIC int postRefineFE(int icalc, NOsh *nosh, FEMparm *feparm, 
  Vfetk *fetk[NOSH_MAXCALC]) {

    int nverts, marked;

    nverts = Gem_numVV(fetk[icalc]->gm);
    if (nverts > feparm->maxvert) {
        Vnm_tprint(1, "    Current number of vertices (%d) exceeds max (%d)!\n",
          nverts, feparm->maxvert);
        return 0;
    }
    Vnm_tprint(1, "      Mesh currently has %d vertices\n", nverts);

    switch(feparm->akeySOLVE) {
        case FRT_UNIF:
            Vnm_tprint(1, "      Commencing uniform refinement.\n");
            break;
        case FRT_GEOM:
            Vnm_tprint(1, "      Commencing geometry-based refinement.\n");
            break;
        case FRT_RESI:
            Vnm_tprint(1, "      Commencing residual-based refinement.\n");
            break;
        case FRT_DUAL:
            Vnm_tprint(1, "      Commencing dual problem-based refinement.\n");
            break;
        case FRT_LOCA:
            Vnm_tprint(1, "      Commencing local-based refinement.\n.");
            break;
        default:
            Vnm_tprint(2, "      Error -- unknown refinement type (%d)!\n", 
              feparm->akeySOLVE);
            return 0;
            break;
    }

    Vnm_print(1, "DEBUG -- akeySOLVE = %d\n", feparm->akeySOLVE);
    marked = AM_markRefine(fetk[icalc]->am, feparm->akeySOLVE, -1, 
      feparm->ekey, feparm->etol);
    if (marked == 0) {
        Vnm_tprint(1, "      Marked 0 simps; hit error/size tolerance.\n");
        return 0;
    }
    Vnm_tprint(1, "      Have %d verts, marked %d.  Refining...\n", nverts,
      marked);
    AM_refine(fetk[icalc]->am, 0, USEHB);
    nverts = Gem_numVV(fetk[icalc]->gm);
    Vnm_tprint(1, "      Done refining; have %d verts.\n", nverts);
    Vnm_redirect(0);
    Gem_shapeChk(fetk[icalc]->gm);
    Vnm_redirect(1);

    return 1;
}


VPUBLIC int writedataFE(int rank, NOsh *nosh, PBEparm *pbeparm, Vfetk *fetk) {

    char writestem[VMAX_ARGLEN];
    char outpath[VMAX_ARGLEN];
    int i, nx, ny, nz, writeit;
    double hx, hy, hzed, xcent, ycent, zcent, xmin, ymin, zmin;
    AM *am;
    Bvec *vec;

    if (nosh->bogus) return 1;

    am = fetk->am;
    vec = am->w0;
  
    for (i=0; i<pbeparm->numwrite; i++) { 

        switch (pbeparm->writetype[i]) {

            writeit = 1;

            case VDT_CHARGE:
 
                Vnm_tprint(2, "    Sorry; can't write charge distribution for FEM!\n");
                writeit = 0;
                break;

            case VDT_POT:
 
                Vnm_tprint(1, "    Writing potential to ");
                Vfetk_fillArray(fetk, vec, VDT_POT);
                break;

            case VDT_SMOL:

                Vnm_tprint(1, "    Writing molecular accessibility to ");
                Vfetk_fillArray(fetk, vec, VDT_SMOL);
                break;

            case VDT_SSPL:

                Vnm_tprint(1, "    Writing spline-based accessibility to ");
                Vfetk_fillArray(fetk, vec, VDT_SSPL);
                break;

            case VDT_VDW:

                Vnm_tprint(1, "    Writing van der Waals accessibility to ");
                Vfetk_fillArray(fetk, vec, VDT_VDW);
                break;

            case VDT_IVDW:

                Vnm_tprint(1, "    Writing ion accessibility to ");
                Vfetk_fillArray(fetk, vec, VDT_IVDW);
                break;

            case VDT_LAP:

                Vnm_tprint(2, "    Sorry; can't write charge distribution for FEM!\n");
                writeit = 0;
                break;

            case VDT_EDENS:

                Vnm_tprint(2, "    Sorry; can't write energy density for FEM!\n");
                writeit = 0;
                break;

            case VDT_NDENS:

                Vnm_tprint(1, "    Writing number density to ");
                Vfetk_fillArray(fetk, vec, VDT_NDENS);
                break;

            case VDT_QDENS:

                Vnm_tprint(1, "    Writing charge density to ");
                Vfetk_fillArray(fetk, vec, VDT_QDENS);
                break;

            case VDT_DIELX:

                Vnm_tprint(2, "    Sorry; can't write x-shifted dielectric map for FEM!\n");
                writeit = 0;
                break;

            case VDT_DIELY:

                Vnm_tprint(2, "    Sorry; can't write y-shifted dielectric map for FEM!\n");
                writeit = 0;
                break;

            case VDT_DIELZ:

                Vnm_tprint(2, "    Sorry; can't write z-shifted dielectric map for FEM!\n");
                writeit = 0;
                break;

            case VDT_KAPPA:

                Vnm_tprint(1, "    Sorry; can't write kappa map for FEM!\n");
                writeit = 0;
                break;

            default:

                Vnm_tprint(2, "Invalid data type for writing!\n");
                writeit = 0;
                return 0;
        }

        if (!writeit) return 0;


#ifdef HAVE_MPI_H
        sprintf(writestem, "%s-PE%d", pbeparm->writestem[i], rank);
#else
        sprintf(writestem, "%s", pbeparm->writestem[i]);
#endif

        switch (pbeparm->writefmt[i]) {

            case VDF_DX:
                sprintf(outpath, "%s.%s", writestem, "dx");
                Vnm_tprint(1, "%s\n", outpath);
                Vfetk_write(fetk, "FILE", "ASC", VNULL, outpath, vec, VDF_DX);
                break;

            case VDF_AVS:
                sprintf(outpath, "%s.%s", writestem, "ucd");
                Vnm_tprint(1, "%s\n", outpath);
                Vfetk_write(fetk, "FILE", "ASC", VNULL, outpath, vec, VDF_AVS);
                break;

            case VDF_UHBD:
                Vnm_tprint(2, "UHBD format not supported for FEM!\n");
                break;

            default:
                Vnm_tprint(2, "Bogus data format (%d)!\n", 
                  pbeparm->writefmt[i]);
                break;
        }
                
    }

    return 1;
}
#endif /* ifdef HAVE_MCX_H */
