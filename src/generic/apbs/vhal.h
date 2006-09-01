/** @defgroup Vhal Vhal class
 *  @brief    A "class" which consists solely of macro definitions which are
 *            used by several other classes
 */

/**
 *  @file       vhal.h
 *  @ingroup    Vhal
 *  @brief      Contains generic macro definitions for APBS
 *  @version  $Id$
 *  @author   Nathan A. Baker
 *
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
 * Copyright (c) 2002-2006.  Washington University in St. Louis.
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
 * Linking APBS statically or dynamically with other modules is making a
 * combined work based on APBS. Thus, the terms and conditions of the GNU
 * General Public License cover the whole combination.
 * 
 * SPECIAL GPL EXCEPTION
 * In addition, as a special exception, the copyright holders of APBS
 * give you permission to combine the APBS program with free software
 * programs and libraries that are released under the GNU LGPL or with
 * code included in releases of ISIM, Ion Simulator Interface, PMV, PyMOL
 * SMOL, VMD, and Vision. Such combined software may be linked with APBS and 
 * redistributed together in original or modified form as mere aggregation
 * without requirement that the entire work be under the scope of the GNU 
 * General Public License. This special exception permission is also extended
 * to any software listed in the SPECIAL GPL EXCEPTION clauses by the PMG,
 * FEtk, MC, or MALOC libraries.
 * 
 * Note that people who make modified versions of APBS are not obligated
 * to grant this special exception for their modified versions; it is
 * their choice whether to do so. The GNU General Public License gives
 * permission to release a modified version without this exception; this
 * exception also makes it possible to release a modified version which
 * carries forward this exception.
 *
 * @endverbatim
 */


#ifndef _VAPBSHAL_H_
#define _VAPBSHAL_H_

/**
 *  @ingroup Vhal
 *  @author  Nathan Baker
 *  @brief   Types of molecular surface definitions
 */
enum eVsurf_Meth {
	VSM_MOL=0, /**<  Ion accessibility is defined using inflated van der Waals
				*    radii, the dielectric coefficient ( ) is defined using the
				*    molecular (Conolly) surface definition without 
                *    smoothing */
    VSM_MOLSMOOTH=1, /**<  As VSM_MOL but with a simple harmonic average
                      *    smoothing */
	VSM_SPLINE=2,    /**<  Spline-based surface definitions. This is primarily
					  * for use with force calculations, since it requires
					  * substantial reparameterization of radii. This is based
					  * on the work of Im et al, Comp. Phys.  Comm. 111 ,
					  * (1998) and uses a cubic spline to define a smoothly
					  * varying characteristic function for the surface-based
					  * parameters. Ion accessibility is defined using inflated
					  * van der Waals radii with the spline function and the
					  * dielectric coefficient is defined using the standard
					  * van der Waals radii with the spline function.  */
	VSM_SPLINE3=3,  /**<  A 5th order polynomial spline is used to create 
					  *  a smoothly varying characteristic function 
					  *  (continuity through 2nd derivatives) for surface
					  *  based paramters. */
	VSM_SPLINE4=4   /**<  A 7th order polynomial spline is used to create 
					  *  a smoothly varying characteristic function 
					  *  (continuity through 3rd derivatives) for surface
					  *  based paramters. */
};

/** @typedef Vsurf_Meth
 *  @ingroup Vhal
 *  @brief   Declaration of the Vsurf_Meth type as the Vsurf_Meth enum
 */
typedef enum eVsurf_Meth Vsurf_Meth;

/**
 * @brief  Version of PBE to solve
 * @ingroup  Vhal
 */
enum eVhal_PBEType {
    PBE_LPBE,  /**<  Traditional Poisson-Boltzmann equation, linearized */
    PBE_NPBE,  /**<  Traditional Poisson-Boltzmann equation, full */
    PBE_LRPBE,  /**<  Regularized Poisson-Boltzmann equation, linearized */
    PBE_NRPBE  /** <  Regularized Poisson-Boltzmann equation, full */
};

/** 
 *  @brief   Declaration of the Vhal_PBEType type as the Vhal_PBEType enum
 *  @ingroup Vhal
 */
typedef enum eVhal_PBEType Vhal_PBEType;


/**
 * @brief Output file format
 * @ingroup Vhal
 */
enum eVoutput_Format {
    OUTPUT_NULL,   /**< No output */
    OUTPUT_FLAT, /**< Output in flat-file format */
    OUTPUT_XML   /**< Output in XML format */
};

/**
 * @brief Declaration of the Voutput_Format type as the VOutput_Format enum
 * @ingroup Vhal
 */
typedef enum eVoutput_Format Voutput_Format;

/**
 * @ingroup  Vhal
 * @author  Nathan Baker
 * @brief  Types of boundary conditions 
 */
enum eVbcfl {
    BCFL_ZERO=0,  /**< Zero Dirichlet boundary conditions */
    BCFL_SDH=1,  /**< Single-sphere Debye-Huckel Dirichlet boundary 
                  * condition */
    BCFL_MDH=2,  /**< Multiple-sphere Debye-Huckel Dirichlet boundary 
                  * condition */
    BCFL_UNUSED=3,  /**< Unused boundary condition method (placeholder) */
    BCFL_FOCUS=4  /**< Focusing Dirichlet boundary condition */
};

/**
 * @brief  Declare Vbcfl type
 * @ingroup  Vhal
 */
typedef enum eVbcfl Vbcfl;

/**
 *  @ingroup Vhal
 *  @author  Nathan Baker
 *  @brief   Types of charge discretization methods
 */
enum eVchrg_Meth {
	VCM_TRIL=0,  /**< Trilinear interpolation of charge to 8 nearest grid
                  *   points.  The traditional method; not particularly good to
                  *   use with PBE forces. */
    VCM_BSPL2=1,  /**< Cubic B-spline across nearest- and
				  *   next-nearest-neighbors.  Mainly for use in grid-sensitive
				  *   applications (such as force calculations). */
    VCM_BSPL4=2  /**< 5th order B-spline for AMOEBA permanent multipoles. */
};

/** @typedef Vchrg_Meth
 *  @ingroup Vhal
 *  @brief   Declaration of the Vchrg_Meth type as the Vchrg_Meth enum
 */
typedef enum eVchrg_Meth Vchrg_Meth;

/**
 *  @ingroup Vhal
 *  @author  Michael Schnieders
 *  @brief   Charge source
 */
enum eVchrg_Src {
	VCM_CHARGE=0,     /**< Partial Charge source distribution */
	VCM_PERMANENT=1,  /**< Permanent Multipole source distribution */
    VCM_INDUCED=2,    /**< Induced Dipole source distribution */
    VCM_NLINDUCED=3   /**< NL Induced Dipole source distribution */
};

/** @typedef Vchrg_Src
 *  @ingroup Vhal
 *  @brief   Declaration of the Vchrg_Src type as the Vchrg_Meth enum
 */
typedef enum eVchrg_Src Vchrg_Src;

/**
 *  @ingroup Vhal
 *  @author  Nathan Baker
 *  @brief   Types of (scalar) data that can be written out of APBS
 */
enum eVdata_Type {
    VDT_CHARGE, /**< Charge distribution (e) */
    VDT_POT,    /**< Potential (kT/e) */
    VDT_SMOL,   /**< Solvent accessibility defined by molecular/Connolly
                 * surface definition (1 = accessible, 0 = inaccessible) */
    VDT_SSPL,   /**< Spline-based solvent accessibility (1 = accessible, 0 =
                 * inaccessible) */
    VDT_VDW,    /**< van der Waals-based accessibility (1 = accessible, 0 =
                 * inaccessible) */
    VDT_IVDW,   /**< Ion accessibility/inflated van der Waals (1 =
                 * accessible, 0 = inaccessible) */
    VDT_LAP,    /**< Laplacian of potential (kT/e/A^2) */
    VDT_EDENS,  /**< Energy density \f$\epsilon (\nabla u)^2\f$, where \f$u\f$
                 * is potential (kT/e/A)^2 */
    VDT_NDENS,  /**< Ion number density \f$\sum c_i \exp (-q_i u)^2\f$, 
	         * where \f$u\f$ is potential (output in M) */
    VDT_QDENS,  /**< Ion charge density \f$\sum q_i c_i \exp (-q_i u)^2\f$,
	         * where \f$u\f$ is potential (output in \f$e_c M\f$) */
    VDT_DIELX,  /**< Dielectric x-shifted map as calculated with the currently
                 * specified scheme (dimensionless) */
    VDT_DIELY,  /**< Dielectric y-shifted map as calculated with the currently
                 * specified scheme (dimensionless) */
    VDT_DIELZ,  /**< Dielectric y-shifted map as calculated with the currently
                 * specified scheme (dimensionless) */
    VDT_KAPPA   /**< Kappa map as calculated with the currently
                 * specified scheme (\f$\AA^{-3}\f$) */
};

/** @typedef Vdata_Type
 *  @ingroup Vhal
 *  @brief   Declaration of the Vdata_Type type as the Vdata_Type enum
 */
typedef enum eVdata_Type Vdata_Type;

/**
 *  @ingroup Vhal
 *  @author  Nathan Baker
 *  @brief   Format of data that can be written out of APBS
 */
enum eVdata_Format {
    VDF_DX=0,  /**< OpenDX (Data Explorer) format */
    VDF_UHBD=1, /**< UHBD format */
    VDF_AVS=2  /**< AVS UCD format */
};

/** @typedef Vdata_Format
 *  @ingroup Vhal
 *  @brief   Declaration of the Vdata_Format type as the Vdata_Format enum
 */
typedef enum eVdata_Format Vdata_Format;

/** 
 * @brief  APBS total execution timer ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_WALL_CLOCK 26

/** 
 * @brief  APBS setup timer ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_SETUP 27

/** 
 * @brief  APBS solver timer ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_SOLVER 28

/** 
 * @brief  APBS energy timer ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_ENERGY 29

/** 
 * @brief  APBS force timer ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_FORCE 30

/** 
 * @brief  APBS temp timer #1 ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_TEMP1 31

/** 
 * @brief  APBS temp timer #2 ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_TEMP2 32

/** @brief The maximum number of molecules that can be involved in a single 
 *         PBE calculation
 *  @ingroup Vhal 
 */
#define MAXMOL 5

/** @brief The maximum number of ion species that can be involved in a single 
 *         PBE calculation
 *  @ingroup Vhal
 */
#define MAXION 10

/** @brief The maximum number of times an MG calculation can be focused
 *  @ingroup Vhal
 */
#define MAXFOCUS 5

/** @brief   Minimum number of levels in a multigrid calculations
 *  @ingroup Vhal
 */
#define VMGNLEV 4

/** @brief   Maximum reduction of grid spacing during a focusing calculation 
 *  @ingroup Vhal
 */ 
#define VREDFRAC 0.25

/** @brief  Number of vertices per simplex (hard-coded to 3D)
 *  @ingroup Vhal
 */
#define VAPBS_NVS 4

/** @brief  Our dimension
 * @ingroup Vhal
 */
#define VAPBS_DIM 3

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH 
 *  @ingroup Vhal
 */
#define VAPBS_RIGHT 0

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH 
 *  @ingroup Vhal
 */
#define VAPBS_FRONT 1

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH 
 *  @ingroup Vhal
 */
#define VAPBS_UP    2

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH 
 *  @ingroup Vhal
 */
#define VAPBS_LEFT  3

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH 
 *  @ingroup Vhal
 */
#define VAPBS_BACK  4

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH 
 *  @ingroup Vhal
 */
#define VAPBS_DOWN  5

/** @brief   A small number used in Vpmg to decide if points are on/off
 *           grid-lines or non-zer0 (etc.)
 *  @ingroup Vhal
 */
#define VPMGSMALL 1e-12


#if defined(VDEBUG)
#   if !defined(APBS_NOINLINE)
#       define APBS_NOINLINE 1
#   endif
#endif

#if !defined(APBS_NOINLINE)

/** @brief   Turns on inlining macros in Vacc class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VACC

/** @brief   Turns on inlining macros in Vatom class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VATOM

/** @brief   Turns on inlining macros in Vcsm class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VCSM

/** @brief   Turns on inlining macros in Vpbe class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VPBE

/** @brief   Turns on inlining macros in Vpee class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VPEE

/** @brief   Turns on inlining macros in Vgreen class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VGREEN

/** @brief   Turns on inlining macros in Vfetk class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VFETK

/** @brief   Turns on inlining macros in Vpmg class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VPMG 

/** @brief   Defines the maximum hash table size in any direction
 *  @ingroup Vhal
 */

#   define MAX_HASH_DIM 75

#endif

/* Fortran name mangling */
#if defined(VF77_UPPERCASE)
#   if defined(VF77_NOUNDERSCORE)
#       define VF77_MANGLE(name,NAME) NAME
#   elif defined(VF77_ONEUNDERSCORE)
#       define VF77_MANGLE(name,NAME) NAME ## _
#   else
#       define VF77_MANGLE(name,NAME) name
#   endif
#else
#   if defined(VF77_NOUNDERSCORE)
#       define VF77_MANGLE(name,NAME) name
#   elif defined(VF77_ONEUNDERSCORE)
#       define VF77_MANGLE(name,NAME) name ## _
#   else
        /** @brief  Name-mangling macro for using FORTRAN functions in C code
         *  @ingroup  Vhal
         */
#       define VF77_MANGLE(name,NAME) name
#   endif
#endif

/* Floating Point Error */
#if defined(MACHINE_EPS)
#     define VFLOOR(value) \
              ((floor(value) != floor(value + MACHINE_EPS)) ? \
              floor(value + MACHINE_EPS) : floor(value))
#else
      /** @brief  Wrapped floor to fix floating point issues in the Intel
       * compiler
       *  @author Todd Dolinksy
       *  @ingroup Vhal
       */
#     define VFLOOR(value) floor(value)
#endif

/* String embedding for ident */
#if defined(HAVE_EMBED)
/**
 * @brief  Allows embedding of RCS ID tags in object files.
 * @author Mike Holst
 * @ingroup Vhal */
#    define VEMBED(rctag) \
         VPRIVATE const char* rctag; \
         static void* use_rcsid=(0 ? &use_rcsid : (void**)&rcsid);
#else
/**
 * @brief  Allows embedding of RCS ID tags in object files.
 * @author Mike Holst
 * @ingroup Vhal */
#    define VEMBED(rctag)
#endif /* if defined(HAVE_EMBED) */

#endif /* #ifndef _VAPBSHAL_H_ */
