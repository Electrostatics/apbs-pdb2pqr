/** @defgroup PBEparm PBEparm class
 *  @brief    Parameter structure for PBE variables independent of solver
 */

/**
 *  @file     pbeparm.h
 *  @ingroup  PBEparm
 *  @brief    Contains declarations for class PBEparm
 *  @version  $Id$
 *  @author   Nathan A. Baker
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
 * Copyright (c) 2002-2003.  Washington University in St. Louis.
 * All Rights Reserved.
 *
 * Portions Copyright (c) 1999-2002.  The Regents of the University of
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

 * @endverbatim
 */

#ifndef _PBEPARM_H_
#define _PBEPARM_H_

/* Generic headers */
#include "maloc/maloc.h"

/* Headers specific to this file */
#include "apbs/vhal.h"

/** @brief   Number of things that can be written out in a single calculation
 *  @ingroup PBEparm
 */
#define PBEPARM_MAXWRITE 10

/**
 *  @struct  PBEparm
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @brief   Parameter structure for PBE variables from input files
 *  @note    If you add/delete/change something in this class, the member
 *           functions -- especially PBEparm_copy -- must be modified
 *           accordingly
 */
struct PBEparm {

    int molid;                 /**< Molecule ID to perform calculation on */
    int setmolid;              /**< Flag, @see molid */
    int useDielMap;            /**< Indicates whether we use an external
                                * dielectric maps (note plural) */
    int dielMapID;             /**< Dielectric map ID (if used) */
    int useKappaMap;           /**< Indicates whether we use an external
                                * kappa map */
    int kappaMapID;            /**< Kappa map ID (if used) */
    int useChargeMap;          /**< Indicates whether we use an external
                                * charge distribution map */
    int chargeMapID;           /**< Charge distribution map ID (if used) */
    int nonlin;                /**< 0 => LPBE, 1 => NPBE */
    int setnonlin;             /**< Flag, @see nonlin */
    int bcfl;                  /**< Boundary condition: 0 => zero, 1 => single
                                * Debye-Huckel sphere, 2 => multiple Debye-
                                * Huckel spheres, 4 => focusing */
    int setbcfl;               /**< Flag, @see bcfl */
    int nion;                  /**< Number of counterion species */
    int setnion;               /**< Flag, @see nion */
    double ionq[MAXION];       /**< Counterion charges (in e) */
    double ionc[MAXION];       /**< Counterion concentrations (in M) */
    double ionr[MAXION];       /**< Counterion radii (in A) */
    int setion[MAXION];        /**< Flag, @see ionq */
    double pdie;               /**< Solute dielectric */
    int setpdie;               /**< Flag, @see pdie */
    double sdie;               /**< Solvent dielectric */
    int setsdie;               /**< Flag, @see sdie */
    Vchrg_Meth chgm;           /**< Charge discretization method */
    int setchgm;               /**< Flag, @see chgm */
    Vsurf_Meth srfm;           /**< Surface calculation method */
    int setsrfm;               /**< Flag, @see srfm */
    double srad;               /**< Solvent radius */
    int setsrad;               /**< Flag, @see srad */
    double swin;               /**< Cubic spline window */
    int setswin;               /**< Flag, @see swin */
    double temp;               /**< Temperature (in K) */
    int settemp;               /**< Flag, @see temp */
    double gamma;              /**< Surface tension for apolar energies/forces
                                * (in kJ/mol/A^2) */
    int setgamma;              /**< Flag, @see gamma */
    int calcenergy;            /**< Energy calculation
				* \li 0: don't calculate out energy
                                * \li 1: calculate total energy 
				* \li 2: calculate atom-level total energies
				*        and all energy components*/
    int setcalcenergy;         /**< Flag, @see calcenergy */
    int calcforce;             /**< Atomic forces I/O 
                                * \li 0: don't calculate forces
                                * \li 1: calculate net forces on molecule
                                * \li 2: calculate atom-level forces */
    int setcalcforce;          /**< Flag, @see calcforce */
    int numwrite;              /**< Number of write statements encountered */
    char writestem[PBEPARM_MAXWRITE][VMAX_ARGLEN]; /**< File stem to write 
                                * data to */
    Vdata_Type writetype[PBEPARM_MAXWRITE]; /**< What data to write */
    Vdata_Format writefmt[PBEPARM_MAXWRITE]; /**< File format to write data 
                                * in */
    int writemat;              /**< Write out the operator matrix? 
                                    \li 0 => no 
                                    \li 1 => yes */
    int setwritemat;           /**< Flag, @see writemat */
    char writematstem[VMAX_ARGLEN];    /**< File stem to write mat */
    int writematflag;          /**< What matrix should we write:
                                    \li 0 => Poisson (differential operator)
                                    \li 1 => Poisson-Boltzmann operator
					     linearized around solution (if
                                             applicable) */


    int parsed;                /**< Has this been filled with anything other
				* than the default values? */
};

/** @typedef PBEparm
 *  @ingroup PBEparm
 *  @brief   Declaration of the PBEparm class as the PBEparm structure
 */
typedef struct PBEparm PBEparm;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (mcsh.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Get charge (e) of specified ion species
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @param   thee PBEparm object
 *  @param   iion Ion species ID
 *  @returns Charge of ion species (e)
 */
VEXTERNC double PBEparm_getIonCharge(PBEparm *thee, int iion);

/** @brief   Get concentration (M) of specified ion species
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @param   thee PBEparm object
 *  @param   iion Ion species ID
 *  @returns Concentration of ion species (M)
 */
VEXTERNC double PBEparm_getIonConc(PBEparm *thee, int iion);

/** @brief   Get radius (A) of specified ion species
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @param   thee PBEparm object
 *  @param   iion Ion species ID
 *  @returns Radius of ion species (A)
 */
VEXTERNC double PBEparm_getIonRadius(PBEparm *thee, int iion);


/** @brief   Construct PBEparm object
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns Newly allocated and initialized PBEparm object
 */
VEXTERNC PBEparm* PBEparm_ctor();

/** @brief   FORTRAN stub to construct PBEparm object
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @param   thee Space for PBEparm object
 *  @returns 1 if succesful, 0 otherwise
 */
VEXTERNC int      PBEparm_ctor2(PBEparm *thee);

/** @brief   Object destructor
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location of PBEparm object
 */
VEXTERNC void     PBEparm_dtor(PBEparm **thee);

/** @brief   FORTRAN stub for object destructor
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to PBEparm object
 */
VEXTERNC void     PBEparm_dtor2(PBEparm *thee);

/** @brief   Consistency check for parameter values stored in object
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @param   thee   PBEparm object
 *  @returns 1 if OK, 0 otherwise
 */
VEXTERNC int      PBEparm_check(PBEparm *thee);

/** @brief   Copy PBEparm object into thee
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @param   thee   PBEparm object (target for copy)
 *  @param   parm   PBEparm object (source for copy)
 */
VEXTERNC void     PBEparm_copy(PBEparm *thee, PBEparm *parm);

/** @brief   Parse a keyword from an input file
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @param   thee   PBEparm object
 *  @param   tok    Token to parse
 *  @param   sock   Stream for more tokens
 *  @return   1 if matched and assigned; -1 if matched, but there's some sort
 *            of error (i.e., too few args); 0 if not matched
 */
VEXTERNC int      PBEparm_parseToken(PBEparm *thee, char tok[VMAX_BUFSIZE],
                    Vio *sock);


#endif 

