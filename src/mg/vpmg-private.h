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
// Internal routines
/////////////////////////////////////////////////////////////////////////// */

/** 
 * @brief  Evaluate a cubic B-spline 
 * @author  Nathan Baker
 * @return  Cubic B-spline value
 */
VPRIVATE double bspline2(
        double x  /** Position */
        );

/** 
 * @brief  Evaluate a cubic B-spline derivative 
 * @author  Nathan Baker
 * @return  Cubic B-spline derivative
 */
VPRIVATE double dbspline2(
        double x  /** Position */
        );

/** 
 * @brief  Calculates charge-potential energy using summation over delta
 *         function positions (i.e. something like an Linf norm)
 * @author  Nathan Baker
 * @return  Energy in kT
 */
VPRIVATE double Vpmg_qfEnergyPoint(
        Vpmg *thee, 
        int extFlag  /** If 1, add external energy contributions to 
                       result */
        );

/** 
 * @brief  Calculates charge-potential energy as integral over a volume
 * @author  Nathan Baker
 * @return  Energy in kT
 */
VPRIVATE double Vpmg_qfEnergyVolume(
        Vpmg *thee, 
        int extFlag  /** If 1, add external energy contributions to 
                       result */
        );

/**
 * @brief  For focusing, fill in the boundaries of the new mesh based on the
 * potential values in the old mesh
 * @author  Nathan Baker
 */
VPRIVATE void focusFillBound(
        Vpmg *thee,  /** New PMG object (the one just created) */
        Vpmg *pmg  /** Old PMG object */
        );

/**
 * @brief  Increment all boundary points by
 *         pre1*(charge/d)*(exp(-xkappa*(d-size))/(1+xkappa*size) to add the
 *         effect of the Debye-Huckel potential due to a single charge
 * @author  Nathan Baker
 */
VPRIVATE void bcfl1(
        double size,  /** Size of the ion */
        double *apos,  /** Position of the ion */
        double charge,  /** Charge of the ion */
        double xkappa,  /** Exponential screening factor */
        double pre1,  /** Unit- and dielectric-dependent prefactor */
        double *gxcf,  /** Set to x-boundary values */
        double *gycf,  /** Set to y-boundary values */
        double *gzcf,  /** Set to z-boundary values */
        double *xf,  /** Boundary point x-coordinates */
        double *yf,  /** Boundary point y-coordinates */
        double *zf,  /** Boundary point z-coordinates */
        int nx,  /** Number of grid points in x-direction */
        int ny,  /** Number of grid points in y-direction */
        int nz /** Number of grid points in y-direction */
        );
/**
 * @brief  Calculate
 *         pre1*(charge/d)*(exp(-xkappa*(d-size))/(1+xkappa*size) due to a
 *         specific ion at a specific point
 * @author  Nathan Baker
 * @returns  Value of above function in arbitrary units (dependent on
 *           pre-factor)
 */
VPRIVATE double bcfl1sp(
        double size,  /** Atom size */
        double *apos,  /** Atom position */
        double charge,  /** Atom charge */ 
        double xkappa,  /** Exponential screening factor */
        double pre1,  /** Unit- and dielectric-dependent prefactor */  
        double *pos  /** Function evaluation position */
        );

/** 
 * @brief  Fill boundary condition arrays
 * @author  Nathan Baker
 */
VPRIVATE void bcCalc(
        Vpmg *thee
        );

/**
 * @brief  Top-level driver to fill all operator coefficient arrays 
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoef(
        Vpmg *thee
        );

/** 
 * @brief  Fill operator coefficient arrays from pre-calculated maps
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefMap(
        Vpmg *thee
        );

/** 
 * @brief  Fill operator coefficient arrays from a molecular surface
 *         calculation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefMol(
        Vpmg *thee
        );

/** 
 * @brief  Fill operator coefficient arrays from a spline-based surface
 *         calculation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefSpline(
        Vpmg *thee
        );

/**
 * @brief  Top-level driver to fill source term charge array
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCharge(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array from a pre-calculated map
 * @author  Nathan Baker
 */
VPRIVATE void fillcoChargeMap(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array from linear interpolation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoChargeSpline1(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array from cubic spline interpolation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoChargeSpline2(
        Vpmg *thee
        );

/**
 * @brief  For focusing, set external energy data members in new Vpmg object
 *         based on energy calculations on old Vpmg object from regions
 *         outside the indicated partition.
 * @author  Nathan Baker, Todd Dolinsky
 */
VPRIVATE void extEnergy(
        Vpmg *thee,  /** Newly created PMG manager */
        Vpmg *pmgOLD,  /** Old PMG manager, source of energies */ 
        PBEparm_calcEnergy extFlag,  /** Energy calculation flag */
        double partMin[3],  /** Partition lower corner */
        double partMax[3],  /** Partition upper corner */ 
        int bflags[6]  /** Which boundaries to include the calculation */
        );

/** 
 * @brief  Charge-field force due to a linear spline charge function
 * @author  Nathan Baker
 */
VPRIVATE void qfForceSpline1(
        Vpmg *thee, 
        double *force,  /** Set to force */
        int atomID  /** Valist atom ID */
        );

/** 
 * @brief  Charge-field force due to a cubic spline charge function
 * @author  Nathan Baker
 */
VPRIVATE void qfForceSpline2(
        Vpmg *thee, 
        double *force,  /** Set to force */
        int atomID  /** Valist atom ID */
        );

#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define IJKx(j,k,i) (((i)*(ny)*(nz))+((k)*(ny))+(j))
#define IJKy(i,k,j) (((j)*(nx)*(nz))+((k)*(nx))+(i))
#define IJKz(i,j,k) (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define VFCHI(iint,iflt) (1.5+((double)(iint)-(iflt)))

/* ///////////////////////////////////////////////////////////////////////////
// External FORTRAN routines 
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

#endif
