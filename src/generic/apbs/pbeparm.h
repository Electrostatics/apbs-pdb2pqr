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

#ifndef _PBEPARM_H_
#define _PBEPARM_H_

#include "apbs/apbs.h"
#include "maloc/maloc.h"

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
                                * dielectric map */
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
    int srfm;                  /**< Surface calculation method
				* \li 0: Mol surface for epsilon; inflated VdW
				*   for kappa; no smoothing 
                                * \li 1: As 0 with harmoic average
                                *   smoothing
                                * \li 2: Cubic spline */
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
                                * \li 2: calculate total energy and all energy
                                *      components*/
    int setcalcenergy;         /**< Flag, @see calcenergy */
    int calcforce;             /**< Atomic forces I/O 
                                * \li 0: don't calculate forces
                                * \li 1: calculate net forces on molecule
                                * \li 2: calculate atom-level forces */
    int setcalcforce;          /**< Flag, @see calcforce */
    int writepot;              /**< 0 => no, 1 => yes */
    int setwritepot;           /**< Flag, @see writepot */
    char writepotstem[VMAX_ARGLEN]; /**< File stem to write pot */
    int writepotfmt;           /**< Potential file formats: 0 => dx, 1 => avs, 
                                * 2 => UHBD */
    int writeacc;              /**< 0 => no, 1 => yes */
    int setwriteacc;           /**< Flag, @see writeacc */
    char writeaccstem[VMAX_ARGLEN];    /**< File stem to write pot */
    int writeaccfmt;           /**< Potential file formats: 0 => dx, 1 => avs, 
                                * 2 => UHBD */
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

