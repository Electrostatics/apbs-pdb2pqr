/** @defgroup Vgreen Vgreen class
 *  @brief    Provides capabilities for pointwise evaluation of free space
 *            Green's function for point charges in a uniform dielectric.
 *  @note     Right now, these are very slow methods without any fast multipole
 *            acceleration.
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

/**
 *  @file     vgreen.h
 *  @ingroup  Vgreen
 *  @brief    Contains declarations for class Vgreen
 *  @version  $Id$
 *  @author   Nathan A. Baker
 */

#ifndef _VGREEN_H_
#define _VGREEN_H_

/* Generic headers */
#include "maloc/maloc.h"
#include "apbs/vhal.h"

/* Specific headers */
#include "apbs/vunit.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"


/**
 *  @struct  Vgreen
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vgreen class/module
 */
struct Vgreen { 

  Valist *alist;      /**< Atom (charge) list for Green's function */
  Vmem *vmem;         /**< Memory management object */
  int initFlagCXXFMM; /**< Flag to indicate whether the C++ FMM object has been
                       * initialized 
                       * @deprecated Not used until we can find a stable 3D FMM
                       * routine */

};

/** @typedef Vgreen
 *  @ingroup Vgreen
 *  @brief   Declaration of the Vgreen class as the Vgreen structure
 */
typedef struct Vgreen Vgreen;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vgreen: Inlineable methods (vgreen.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VGREEN)

    /** @brief   Get the atom list associated with this Green's function object
     *  @ingroup Vgreen
     *  @author  Nathan Baker
     *  @param   thee  Vgreen object
     *  @return  Pointer to Valist object associated with this Green's function
     *           object
     */
    VEXTERNC Valist* Vgreen_getValist(Vgreen *thee);

    /** @brief   Return the memory used by this structure (and its contents)
     *           in bytes
     *  @ingroup Vgreen
     *  @author  Nathan Baker
     *  @param   thee  Vgreen object
     *  @return  The memory used by this structure and its contents in bytes
     */
    VEXTERNC int     Vgreen_memChk(Vgreen *thee);

#else /* if defined(VINLINE_VGREEN) */
#   define Vgreen_getValist(thee) ((thee)->alist)
#   define Vgreen_memChk(thee) (Vmem_bytes((thee)->vmem))
#endif /* if !defined(VINLINE_VGREEN) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vgreen: Non-Inlineable methods (vgreen.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct the Green's function oracle
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   alist  Atom (charge) list associated with object
 *  @return  Pointer to newly allocated Green's function oracle
 */
VEXTERNC Vgreen* Vgreen_ctor(Valist *alist);

/** @brief   FORTRAN stub to construct the Green's function oracle
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory allocated for object
 *  @param   alist  Atom (charge) list associated with object
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int     Vgreen_ctor2(Vgreen *thee, Valist *alist);

/** @brief   Destruct the Green's function oracle
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   thee Pointer to memory location for object
 */
VEXTERNC void    Vgreen_dtor(Vgreen **thee);

/** @brief   FORTRAN stub to destruct the Green's function oracle
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   thee Pointer to object
 */
VEXTERNC void    Vgreen_dtor2(Vgreen *thee);

/** @brief   Get the Green's function for Helmholtz's equation integrated over
 *           the atomic point charges
 *  
 *           Returns the potential \f$\phi\f$ defined by
 *           \f[ \phi(r) = \sum_i \frac{q_i e^{-\kappa r_i}}{r_i} \f]
 * 
 *           where \f$\kappa\f$ is the inverse screening length (in \f$\AA\f$),
 *           \f$q_i\f$ is the atomic charge (in e), and \f$r_i\f$ r_i is the
 *           distance from atom \f$i\f$ to the observation point \f$r\f$.  The
 *           potential is scaled to units of V.
 *
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @bug     Not implemented yet
 *  @note    Not implemented yet
 *  @param   thee  Vgreen object
 *  @param   position  The coordinates of \f$r\f$ (see above)
 *  @param   dim  The dimension of the space in which this is evaluated
 *  @param   kappa The value of \f$\kappa\f$ (see above)
 *  @return  The potential value as defined above in units of V
 */
VEXTERNC double  Vgreen_helmholtz(Vgreen *thee, double *position, double dim,
                   double kappa);

/** @brief   Get the gradient of Green's function for Helmholtz's equation
 *           integrated over the atomic point charges
 *
 *           Returns the field \f$\nabla \phi\f$ defined by
 *           \f[ \nabla \phi(r) = \nabla \sum_i \frac{q_i e^{-\kappa r_i}}{r_i}
 *           \f]
 *
 *           where \f$\kappa\f$ is the inverse screening length (in \f$\AA\f$),
 *           \f$q_i\f$ is the atomic charge (in e), and \f$r_i\f$ r_i is the
 *           distance from atom \f$i\f$ to the observation point \f$r\f$.  The
 *           potential is scaled to units of V/\f$\AA\f$.
 *
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @bug     Not implemented yet
 *  @note    Not implemented yet
 *  @param   thee  Vgreen object
 *  @param   position  The coordinates of \f$r\f$ (see above)
 *  @param   dim  The dimension of the space in which this is evaluated
 *  @param   kappa The value of \f$\kappa\f$ (see above)
 *  @param   grad The field as defined above in units of V/\f$\AA\f$
 */
VEXTERNC void    Vgreen_helmholtzD(Vgreen *thee, double *position, 
                   double dim, double kappa, double *grad);

/** @brief   Get the Coulomb's Law Green's function (solution to Laplace's
 *           equation) integrated over the atomic point charges
 * 
 *           Returns the potential \f$\phi\f$ defined by 
 *           \f[ \phi(r) = \sum_i \frac{q_i}{r_i} \f]
 *           where \f$q_i\f$ is the atomic charge (in e) and \f$r_i\f$ is the
 *           distance to the observation point \f$r\f$.  The potential is
 *           scaled to units of V.
 *
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   thee Vgreen object
 *  @param   position  The coordinates of \f$r\f$ (see above)
 *  @param   dim  The dimension of the problem domain
 *  @return  The potential defined above in units of V
 */
VEXTERNC double  Vgreen_coulomb(Vgreen *thee, double *position, double dim);

/** @brief   Get gradient of the Coulomb's Law Green's function (solution to
 *           Laplace's equation) integrated over the atomic point charges
 *
 *           Returns the field \f$\nabla \phi\f$ defined by
 *           \f[ \nabla \phi(r) = \sum_i \frac{q_i}{r_i} \f]
 *           where \f$q_i\f$ is the atomic charge (in e) and \f$r_i\f$ is the
 *           distance to the observation point \f$r\f$.  The field is
 *           scaled to units of V/\f$\AA\f$.
 *
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   thee Vgreen object
 *  @param   position  The coordinates of \f$r\f$ (see above)
 *  @param   dim  The dimension of the problem domain
 *  @param   grad The field defined above in units of V
 */
VEXTERNC void    Vgreen_coulombD(Vgreen *thee, double *position, double dim,
                   double *grad);

#endif /* ifndef _VGREEN_H_ */
