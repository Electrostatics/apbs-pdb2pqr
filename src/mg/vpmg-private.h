/**
 *  @file    vpmg-private.h
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @brief   Class Vpmg private method declaration
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

#ifndef _VPMG_PRIVATE_H_
#define _VPMG_PRIVATE_H_

#include "apbscfg.h"
#include "apbs/vpmg.h"

/* ///////////////////////////////////////////////////////////////////////////
// Private routines
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double bspline2(double x);
VPUBLIC double dbspline2(double x);

/* ///////////////////////////////////////////////////////////////////////////
// External FORTRAN ROUTINES 
/////////////////////////////////////////////////////////////////////////// */
#define F77BCOLCOMP VF77_MANGLE(bcolcomp, BCOLCOMP)
VEXTERNC void F77BCOLCOMP(int *iparm, double *rparm, int *iwork, 
  double *rwork, double *nzval, int *rowind, int *colptr, int *flag);

#define F77PCOLCOMP VF77_MANGLE(pcolcomp, PCOLCOMP)
VEXTERNC void F77PCOLCOMP(int *nrow, int *ncol, int *nonz, 
  double *nzval, int *rowind, int *colptr, 
  char *path, char *title, char *mxtype);

#define F77MGSZ VF77_MANGLE(mgsz, MGSZ)
VEXTERNC void F77MGSZ(int *mgcoar, int *mgdisc, int *mgsolv, int *nx, int *ny,
  int *nz, int *nlev, int *nxc, int *nyc, int *nyz, int *nf, int *nc, 
  int *narr, int *narrc, int *n_rpc, int *n_iz, int *n_ipc, int *nrwk, 
  int *niwk);

#define F77PACKMG VF77_MANGLE(packmg, PACKMG)
VEXTERNC void F77PACKMG(int *iparm, double *rparm, int *nrwk, int *niwk,
  int *nx, int *ny, int *nz, int *nlev, int *nu1, int *nu2, int *mgkey, 
  int *itmax, int *istop, int *ipcon, int *nonlin, int *mgsmoo, int *mgprol, 
  int *mgcoar, int *mgsolv, int *mgdisc, int *iinfo, double *errtol,
  int *ipkey, double *omegal, double *omegan, int *irite, int *iperf);

#define F77CGMGDRIV VF77_MANGLE(cgmgdriv, CGMGDRIV)
VEXTERNC void F77CGMGDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77NEWDRIV    VF77_MANGLE(newdriv, NEWDRIV)
VEXTERNC void F77NEWDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77MGDRIV     VF77_MANGLE(mgdriv, MGDRIV)
VEXTERNC void F77MGDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77NCGHSDRIV  VF77_MANGLE(ncghsdriv, NCGHSDRIV)
VEXTERNC void F77NCGHSDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77NSORDRIV   VF77_MANGLE(nsordriv, NSORDRIV)
VEXTERNC void F77NSORDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77NGSRBDRIV  VF77_MANGLE(ngsrbdriv, NGSRBDRIV)
VEXTERNC void F77NGSRBDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77NWJACDRIV  VF77_MANGLE(nwjacdriv, NWJACDRIV)
VEXTERNC void F77NWJACDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77NRICHDRIV  VF77_MANGLE(nrichdriv, NRICHDRIV)
VEXTERNC void F77NRICHDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77TSECND     VF77_MANGLE(tsecnd, TSECND)
#define F77VPMGANORM  VF77_MANGLE(vpmganorm, VPMGANORM)
#define F77VPMGABAND  VF77_MANGLE(vpmgaband, VPMGABAND)
#define F77DPBFA      VF77_MANGLE(dpbfa, DPBFA)
#define F77DPBDI      VF77_MANGLE(dpbdi, DPBDI)
#define F77EIGDRIV    VF77_MANGLE(eigdriv, EIGDRIV)
#define F77ANORMDRIV  VF77_MANGLE(anormdriv, ANORMDRIV)

#define F77MYPDEFINIT VF77_MANGLE(mypdefinit, MYPDEFINIT)
VEXTERNC void F77MYPDEFINIT(int *nion, double *ionQ, double *ionConc);

#define F77MYPDEFCLEAR VF77_MANGLE(mypdefclear, MYPDEFCLEAR)
VEXTERNC void F77MYPDEFCLEAR();

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmg: Private methods
/////////////////////////////////////////////////////////////////////////// */
#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define IJKx(j,k,i) (((i)*(ny)*(nz))+((k)*(ny))+(j))
#define IJKy(i,k,j) (((j)*(nx)*(nz))+((k)*(nx))+(i))
#define IJKz(i,j,k) (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define VFCHI(iint,iflt) (1.5+((double)(iint)-(iflt)))

#endif
