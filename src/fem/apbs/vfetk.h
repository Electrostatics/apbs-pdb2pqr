/** @defgroup Vfetk Vfetk class
 *  @brief    FEtk master class (interface between FEtk and APBS)
 */

/**
 *  @file     vfetk.h
 *  @ingroup  Vfetk
 *  @brief    Contains declarations for class Vfetk
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
 * Copyright (c) 2002.  Washington University in St. Louis.
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

#ifndef _VFETK_H_
#define _VFETK_H_

#include "maloc/maloc.h"
#include "mc/mc.h"
#include "apbs/vhal.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vcsm.h"
#include "apbs/vpbe.h"
#include "apbs/vunit.h"
#include "apbs/vgreen.h"
#include "apbs/vcap.h"

/**
 *  @struct  Vfetk
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vfetk class/module
 *
 *  Many of the routines and macros are borrowed from the main.c driver
 *  (written by Mike Holst) provided with the PMG code.
 *
 */
struct Vfetk { 

  Vmem *vmem;         /**< Memory management object */
  Gem *gm;            /**< Grid manager (container class for master vertex
                       * and simplex lists as well as prolongation
                       * operator for updating after refinement ) */
  AM *am;             /**< Multilevel algebra manager */

  Vpbe *pbe;          /**< Poisson-Boltzmann object */
  Vcsm *csm;          /**< Charge-simplex map */
  
};

/** @typedef Vfetk
 *  @ingroup Vfetk
 *  @brief   Declaration of the Vfetk class as the Vfetk structure
 */
typedef struct Vfetk Vfetk;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vfetk: Inlineable methods (vfetk.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VFETK)

    /** @brief   Get a pointer to the Gem (grid manager) object
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @param   thee  Vfetk object
     *  @return  Pointer to the Gem (grid manager) object
     */
    VEXTERNC Gem*    Vfetk_getGem(Vfetk *thee);

    /** @brief   Get a pointer to the AM (algebra manager) object
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @param   thee  Vfetk object
     *  @return  Pointer to the AM (algebra manager) object
     */
    VEXTERNC AM*     Vfetk_getAM(Vfetk *thee);

    /** @brief   Get a pointer to the Vpbe (PBE manager) object
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @param   thee  Vfetk object
     *  @return  Pointer to the Vpbe (PBE manager) object
     */
    VEXTERNC Vpbe*   Vfetk_getVpbe(Vfetk *thee);

    /** @brief   Get a pointer to the Vcsm (charge-simplex map) object
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @param   thee  Vfetk object
     *  @return  Pointer to the Vcsm (charge-simplex map) object
     */
    VEXTERNC Vcsm*   Vfetk_getVcsm(Vfetk *thee);

    /** @brief   Get the partition information for a particular atom
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @note    Friend function of Vatom
     *  @param   thee  Vfetk object
     *  @param   iatom Valist atom ID
     *  @returns Partition ID 
     */
    VEXTERNC int     Vfetk_getAtomColor(Vfetk *thee, int iatom);

#else /* if defined(VINLINE_VFETK) */
#   define Vfetk_getGem(thee) ((thee)->gm)
#   define Vfetk_getAM(thee) ((thee)->am)
#   define Vfetk_getVpbe(thee) ((thee)->pbe)
#   define Vfetk_getVcsm(thee) ((thee)->csm)
#   define Vfetk_getAtomColor(thee, iatom) (Vatom_getPartID(Valist_getAtom(Vpbe_getValist(thee->pbe), iatom)))
#endif /* if !defined(VINLINE_VFETK) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vfetk: Non-Inlineable methods (vfetk.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Constructor for Vfetk object
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param   pbe  Vpbe (PBE manager) object
 *  @param   gm   Gem (geometry manager) object
 *  @param   am   AM (algebra manager) object
 *  @return  Pointer to newly allocated Vfetk object 
 */
VEXTERNC Vfetk*  Vfetk_ctor(Vpbe *pbe, Gem *gm, AM *am);

/** @brief   FORTRAN stub constructor for Vfetk object
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param   thee Vfetk obeject memory address
 *  @param   apbe  Vpbe (PBE manager) object
 *  @param   gm   Gem (geometry manager) object
 *  @param   am   AM (algebra manager) object
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int     Vfetk_ctor2(Vfetk *thee, Vpbe *apbe, Gem *gm, AM *am);

/** @brief   Object destructor
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void    Vfetk_dtor(Vfetk **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void    Vfetk_dtor2(Vfetk *thee);

/** @brief   Create an array containing the solution (electrostatic potential
 *           in units of \f$k_B T/e\f$) at the finest mesh level. 
 *  @ingroup Vfetk
 *  @author  Nathan Baker and Michael Holst
 *  @note    The user is responsible for destroying the newly created array
 *  @param   thee   Vfetk object
 *  @param   length Set to the length of the newly created array
 *  @return  Newly created array of length "length" (see above); the user is
 *           responsible for destruction
 */
VEXTERNC double* Vfetk_getSolution(Vfetk *thee, int *length);

/** @brief   Return the total electrostatic energy
 * 
 *   Using the solution at the finest mesh level, get the electrostatic energy
 *   using the free energy functional for the Poisson-Boltzmann equation
 *   without removing any self-interaction terms (i.e., removing the reference
 *   state of isolated charges present in an infinite dielectric continuum with
 *   the same relative permittivity as the interior of the protein) and return
 *   the result in units of \f$k_B T\f$.  The argument color allows the user to
 *   control the partition on which this energy is calculated; if (color == -1)
 *   no restrictions are used.  The solution is obtained from the finest level
 *   of the passed AM object, but atomic data from the Vfetk object is used to
 *   calculate the energy.
 * 
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param thee Vfetk object
 *  @param color Partition restriction; if non-negative, energy calculation is
 *               restricted to the specified partition.
 *  @param nonlin If 1, the NPBE energy functional is used, if 0 then the LPBE
 *                energy functional is used.
 *  @return Total electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double  Vfetk_energy(Vfetk *thee, int color, int nonlin);

/** @brief   Get the "mobile charge" and "polarization" contributions to the
 *           electrostatic energy.
 * 
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the mobile charges
 *           with the potential and polarization of the dielectric medium:
 *              \f[ G = \frac{1}{4 I_s} \sum_i c_i q_i^2 \int
 *              \overline{\kappa}^2(x) e^{-q_i u(x)} dx + \frac{1}{2} \int 
 *              \epsilon ( \nabla u )^2 dx \f]
 *           for the NPBE and
 *              \f[ G = \frac{1}{2} \int \overline{\kappa}^2(x) u^2(x) dx + 
 *              \frac{1}{2} \int \epsilon ( \nabla u )^2 dx \f]
 *           for the LPBE.  Here \f$i\f$ denotes the counterion species,
 *           \f$I_s\f$ is the bulk ionic strength, \f$\overline{\kappa}^2(x)\f$
 *           is the modified Debye-Huckel parameter, \f$c_i\f$ is the  
 *           concentration of species \f$i\f$, \f$q_i\f$ is the charge of
 *           species \f$i\f$, \f$\epsilon\f$ is the dielectric function, and
 *           \f$u(x)\f$ is the dimensionless electrostatic potential.  The
 *           energy is scaled to units of \f$k_b T\f$.
 *
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param thee  Vfetk object
 *  @param color Partition restriction for energy evaluation, only used if
 *               non-negative
 *  @return The "mobile charge" and "polarization" contributions to the
 *           electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double  Vfetk_dqmEnergy(Vfetk *thee, int color);

/** @brief   Get the "fixed charge" contribution to the electrostatic energy
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the fixed charges
 *           with the potential: \f[ G = \sum_i q_i u(r_i) \f]
 *           and return the result in units of \f$k_B T\f$.  Clearly, no
 *           self-interaction terms are removed.  A factor a 1/2 has to be
 *           included to convert this to a real energy.
 *
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param   thee   Vfetk object
 *  @param   color Partition restriction for energy evaluation, only used if
 *               non-negative
 *  @returns The fixed charge electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double  Vfetk_qfEnergy(Vfetk *thee, int color);

/** @brief   Calculate the log determinant of the specified operator.
 *  @ingroup Vfetk
 *  @author  Nathan Baker and Stephen Bond
 *  @note    \li Only works with symmetric positive definite matrices
 *           \li Uses LU or Recycled Cholesky factorization and can be 
 *           very memory- and time-consuming for large matrices.
 *  @param   thee  Vfetk object
 *  @param   color Partition to evaluate over (ignored if <0)
 *  @param   oflag  Operator to evaluate:
 *           \li 0:  Helmholtz operator (NPBE tangent operator evaluated at zero
 *                 solution)
 *           \li 1:  Response function (NPBE tangent operator evaluated at NPBE
 *                 solution)
 *  @param   mflag  Method to use:
 *           \li 0:  Full nonsymmetric SuperLU factor with ROW/COL reordering
 *           \li 1:  Recycled symmetric Cholesky factor with no reordering
 *  @return  The log determinant of the specified operator
 *  @bug     color argument ignored
 */
VEXTERNC double Vfetk_lnDet(Vfetk *thee, int color, int oflag, int mflag);

/** @brief   Return the memory used by this structure (and its contents)
 *           in bytes
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param   thee  Vfetk object
 *  @return  The memory used by this structure and its contents in bytes
 */
VEXTERNC int     Vfetk_memChk(Vfetk *thee);

/** @brief   Transfer color (partition ID) information frmo a partitioned mesh
 *           to the atoms.
 * 
 *           Transfer color information from partitioned mesh to the atoms.
 *           In the case that a charge is shared between two partitions, the
 *           partition color of the first simplex is selected.  Due to the
 *           arbitrary nature of this selection, THIS METHOD SHOULD ONLY BE
 *           USED IMMEDIATELY AFTER PARTITIONING!!!
 *  @warning This function should only be used immediately after mesh
 *           partitioning
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @note    This is a friend function of Vcsm
 *  @param   thee Vfetk object
 */
VEXTERNC void    Vfetk_setAtomColors(Vfetk *thee);

/** @brief   Writes a Bmat to disk in Harwell-Boeing sparse matrix format.
 * 
 *  @ingroup Vfetk
 *  @author  Stephen Bond
 *  @note    This is a friend function of Bmat
 *  @param   thee Bmat object
 *  @param   fname char Output filename
 *  @bug     Hardwired to only handle the single block symmetric case.
 */
VEXTERNC void    Bmat_printHB(Bmat *thee, char *fname);

/** @brief   Assembles the Cholesky factorization of a Bmat.
 * 
 *  @ingroup Vfetk
 *  @author  Stephen Bond
 *  @note    This is a friend function of Bmat
 *  @param   thee Bmat object
 *  @param   flag  Type of factor to be assembled:
 *           \li 0:  Full Cholesky factor
 *           \li 1:  Diagonal of the Cholesky factor
 */
VEXTERNC int     Bmat_choleskyFactor(Bmat *thee, int flag);

/** @brief   Returns the log(abs(det(D))) of a diagonal Mat, D.
 * 
 *  @ingroup Vfetk
 *  @author  Stephen Bond
 *  @note    This is a friend function of Mat
 *  @param   thee Mat object
 *  @bug     Only works for Mat's of type RLN, CLN, or DRC.
 */
VEXTERNC double  Mat_lnDetDiag(Mat *thee);

#endif /* ifndef _VFETK_H_ */
