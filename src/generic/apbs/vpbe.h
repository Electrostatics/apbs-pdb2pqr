/** @defgroup Vpbe Vpbe class
 *  @brief    The Poisson-Boltzmann master class
 *    
 *            Contains objects and parameters used in every PBE calculation,
 *            regardless of method.
 * 
 */

/**
 *  @file       vpbe.h
 *  @ingroup    Vpbe
 *  @brief      Contains declarations for class Vpbe
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
 * Copyright (c) 1999-2002.  The Regents of the University of California.
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

#ifndef _VPBE_H_
#define _VPBE_H_

#include "maloc/maloc.h"
#include "apbs/vhal.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vacc.h"
#include "apbs/vunit.h"
#include "apbs/vgreen.h"

/**
 *  @struct  Vpbe
 *  @ingroup Vpbe
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vpbe class/module
 */
struct Vpbe { 

  Vmem *vmem;         /**< Memory management object */

  Valist *alist;      /**< Atom (charge) list */
  Vacc *acc;          /**< Accessibility object */
  Vgreen *green;      /**< Green's function oracle */

  double T;           /**< Temperature (K) */
  double soluteDiel;  /**< Solute dielectric constant (unitless) */
  double solventDiel; /**< Solvent dielectric constant (unitless) */
  double solventRadius;
                      /**< Solvent probe radius (angstroms) for accessibility;
                       * determining defining volumes for the dielectric
                       * coefficient */
  double bulkIonicStrength; /**< Bulk ionic strength (M) */
  double maxIonRadius;      /**< Max ion radius (A; used for calculating
                             * accessiblity and defining volumes for ionic
                             * strength coeffcients) */
  int    numIon;            /**< Total number of ion species */
  double ionConc[MAXION];   /**< Concentration (M) of each species */
  double ionRadii[MAXION];  /**< Ionic radius (A) of each species */
  double ionQ[MAXION];      /**< Charge (e) of each species */

  double xkappa;      /**< Debye-Huckel parameter (bulk) */
  double deblen;      /**< Debye length (bulk) */
  double zkappa2;     /**< Square of modified Debye-Huckel parameter (bulk) */
  double zmagic;      /**< Delta function scaling parameter */

  double soluteCenter[3]; /**< Center of solute molecule (A) */
  double soluteRadius; /**< Radius of solute molecule (A) */
  double soluteXlen;  /**< Solute length in x-direction */
  double soluteYlen;  /**< Solute length in y-direction */
  double soluteZlen;  /**< Solute length in z-direction */
  double soluteCharge; /**< Charge of solute molecule (e) */

  int paramFlag;      /**< Check to see if the parameters have been set */
  
};

/** @typedef Vpbe
 *  @ingroup Vpbe
 *  @brief   Declaration of the Vpbe class as the Vpbe structure
 */
typedef struct Vpbe Vpbe;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Inlineable methods (vpbe.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPBE)

    /** @brief   Get atom list
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Pointer to internal Valist object
     */
    VEXTERNC Valist* Vpbe_getValist(Vpbe *thee);

    /** @brief   Get accessibility oracle
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Pointer to internal Vacc object
     */
    VEXTERNC Vacc*   Vpbe_getVacc(Vpbe *thee);

    /** @brief   Get Green's function oracle
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Pointer to internal Vgreen object
     */
    VEXTERNC Vgreen* Vpbe_getVgreen(Vpbe *thee);

    /** @brief   Get bulk ionic strength
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Bulk ionic strength (M)
     */
    VEXTERNC double  Vpbe_getBulkIonicStrength(Vpbe *thee);

    /** @brief   Get maximum radius of ion species
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Maximum radius (A)
     */
    VEXTERNC double  Vpbe_getMaxIonRadius(Vpbe *thee);

    /** @brief   Get temperature
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Temperature (K)
     */
    VEXTERNC double  Vpbe_getTemperature(Vpbe *thee);           

    /** @brief   Get solute dielectric constant
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Solute dielectric constant
     */
    VEXTERNC double  Vpbe_getSoluteDiel(Vpbe *thee); 

    /** @brief   Get sphere radius which bounds biomolecule
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Sphere radius which bounds biomolecule (A)
     */
    VEXTERNC double  Vpbe_getSoluteRadius(Vpbe *thee);

    /** @brief   Get length of solute in x dimension
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Length of solute in x dimension (A)
     */
    VEXTERNC double  Vpbe_getSoluteXlen(Vpbe *thee);

    /** @brief   Get length of solute in y dimension
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Length of solute in y dimension (A)
     */
    VEXTERNC double  Vpbe_getSoluteYlen(Vpbe *thee);

    /** @brief   Get length of solute in z dimension
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Length of solute in z dimension (A)
     */
    VEXTERNC double  Vpbe_getSoluteZlen(Vpbe *thee);

    /** @brief   Get coordinates of solute center
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Pointer to 3*double array with solute center coordinates (A)
     */
    VEXTERNC double* Vpbe_getSoluteCenter(Vpbe *thee);

    /** @brief   Get total solute charge
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Total solute charge (e)
     */
    VEXTERNC double  Vpbe_getSoluteCharge(Vpbe *thee);

    /** @brief   Get solvent dielectric constant
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Solvent dielectric constant
     */
    VEXTERNC double  Vpbe_getSolventDiel(Vpbe *thee);

    /** @brief   Get solvent molecule radius
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Solvent molecule radius (A)
     */
    VEXTERNC double  Vpbe_getSolventRadius(Vpbe *thee);

    /** @brief   Get Debye-Huckel parameter 
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Bulk Debye-Huckel parameter (\f$\AA^{-1}\f$)
     */
    VEXTERNC double  Vpbe_getXkappa(Vpbe *thee);

    /** @brief   Get Debye-Huckel screening length
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Debye-Huckel screening length (\f$\AA\f$)
     */
    VEXTERNC double  Vpbe_getDeblen(Vpbe *thee);

    /** @brief   Get modified squared Debye-Huckel parameter
     *  @ingroup Vpbe
     *  @author  Nathan Baker
     *  @param   thee Vpbe object
     *  @return  Modified squared Debye-Huckel parameter (\f$\AA^{-2}\f$)
     */
    VEXTERNC double  Vpbe_getZkappa2(Vpbe *thee);

    /** @brief   Get charge scaling factor
     *  @ingroup Vpbe
     *  @author  Nathan Baker and Mike Holst
     *  @param   thee Vpbe object
     *  @return  Get factor for scaling charges (in e) to internal units
     */
    VEXTERNC double  Vpbe_getZmagic(Vpbe *thee);

#else /* if defined(VINLINE_VPBE) */
#   define Vpbe_getValist(thee) ((thee)->alist)
#   define Vpbe_getVacc(thee) ((thee)->acc)
#   define Vpbe_getVgreen(thee) ((thee)->green)
#   define Vpbe_getBulkIonicStrength(thee) ((thee)->bulkIonicStrength)
#   define Vpbe_getTemperature(thee) ((thee)->T)           
#   define Vpbe_getSoluteDiel(thee) ((thee)->soluteDiel) 
#   define Vpbe_getSoluteCenter(thee) ((thee)->soluteCenter)
#   define Vpbe_getSoluteRadius(thee) ((thee)->soluteRadius)
#   define Vpbe_getSoluteXlen(thee) ((thee)->soluteXlen)
#   define Vpbe_getSoluteYlen(thee) ((thee)->soluteYlen)
#   define Vpbe_getSoluteZlen(thee) ((thee)->soluteZlen)
#   define Vpbe_getSoluteCharge(thee) ((thee)->soluteCharge)
#   define Vpbe_getSolventDiel(thee) ((thee)->solventDiel)
#   define Vpbe_getSolventRadius(thee) ((thee)->solventRadius)
#   define Vpbe_getMaxIonRadius(thee) ((thee)->maxIonRadius)
#   define Vpbe_getXkappa(thee) ((thee)->xkappa)
#   define Vpbe_getDeblen(thee) ((thee)->deblen)
#   define Vpbe_getZkappa2(thee) ((thee)->zkappa2)
#   define Vpbe_getZmagic(thee) ((thee)->zmagic)
#endif /* if !defined(VINLINE_VPBE) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Non-Inlineable methods (vpbe.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct Vpbe object
 *  @ingroup Vpbe
 *  @author  Nathan Baker and Mike Holst
 *  @notes   This is partially based on some of Mike Holst's PMG code.  Here
 *           are a few of the original function comments:
 *           kappa is defined as follows:
 *           \f[ \kappa^2 = \frac{8 \pi N_A e_c^2 I_s}{1000 \epsilon_w k_B T}
 *           \f] 
 *           where the units are esu*esu/erg/mol.  To obtain \f$\AA^{-2}\f$, we
 *           multiply by \f$10^{-16}\f$.  Thus, in \f$\AA^{-2}\f$, where
 *           \f$k_B\f$ and \f$e_c\f$ are in gaussian rather than mks units, the
 *           proper value for kappa is:
 *           \f[ \kappa^2 = \frac{8 \pi N_A e_c^2 I_s}{1000 \epsilon_w k_b T}
 *           \times 10^{-16} \f]
 *           and the factor of \f$10^{-16}\f$ results from converting cm^2 to
 *           angstroms^2, noting that the 1000 in the denominator has converted
 *           m^3 to cm^3, since the ionic strength \f$I_s\f$ is assumed to have
 *           been provided in moles per liter, which is moles per 1000 cm^3.
 *  @param   alist  Atom list
 *  @param   ionNum  Number of counterion species
 *  @param   ionConc Array containing counterion species' concentrations (M)
 *  @param   ionRadii Array containing counterion species' radii (A)
 *  @param   ionQ Array containing counterion species' charges (e)
 *  @param   T temperature (K)
 *  @param   soluteDiel Solute dielectric constant
 *  @param   solventDiel Solvent dielectric constant
 *  @param   solventRadius Solvent radius
 *  @return  Pointer to newly allocated Vpbe object
 */
VEXTERNC Vpbe*   Vpbe_ctor(Valist *alist, int ionNum, double *ionConc, 
		    double *ionRadii, double *ionQ, double T,
                    double soluteDiel, double solventDiel,  
                    double solventRadius);

/** @brief   FORTRAN stub to construct Vpbe objct
 *  @ingroup Vpbe
 *  @author  Nathan Baker and Mike Holst
 *  @notes   This is partially based on some of Mike Holst's PMG code.  Here
 *           are a few of the original function comments:
 *           kappa is defined as follows:
 *           \f[ \kappa^2 = \frac{8 \pi N_A e_c^2 I_s}{1000 eps_w k_B T} \f]
 *           where the units are esu*esu/erg/mol.  To obtain \f$\AA^{-2}\f$, we
 *           multiply by \f$10^{-16}\f$.
 *           Thus, in \f$\AA^{-2}\f$, where \f$k_B\f$ and \f$e_c\f$ are in 
 *           gaussian rather than mks units, the proper value for kappa is:
 *           \f[ \kappa^2 = \frac{8 pi N_A e_c^2 I_s}{1000 eps_w k_b T} \times 
 *           10^{-16} \f]
 *           and the factor of \f$10^{-16}\f$ results from converting cm^2 to 
 *           angstroms^2, noting that the 1000 in the denominator has converted
 *           m^3 to cm^3, since the ionic strength \f$I_s\f$ is assumed to have
 *           been provided in moles per liter, which is moles per 1000 cm^3. 
 *  @param   thee   Pointer to memory allocated for Vpbe object
 *  @param   alist  Atom list
 *  @param   ionNum  Number of counterion species
 *  @param   ionConc Array containing counterion species' concentrations (M)
 *  @param   ionRadii Array containing counterion species' radii (A)
 *  @param   ionQ Array containing counterion species' charges (e)
 *  @param   T temperature (K)
 *  @param   soluteDiel Solute dielectric constant
 *  @param   solventDiel Solvent dielectric constant
 *  @param   solventRadius Solvent radius
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int    Vpbe_ctor2(Vpbe *thee, Valist *alist, int ionNum, 
		    double *ionConc, double *ionRadii, double *ionQ, 
                    double T, double soluteDiel, 
                    double solventDiel, double solventRadius);

/** @brief   Get information about the counterion species present
 *  @ingroup Vpbe
 *  @author  Nathan Baker
 *  @param   thee   Pointer to Vpbe object
 *  @param   nion   Set to the number of counterion species
 *  @param   ionConc Array to store counterion species' concentrations (M)
 *  @param   ionRadii Array to store counterion species' radii (A)
 *  @param   ionQ Array to store counterion species' charges (e)
 *  @return  Number of ions
 */
VEXTERNC int     Vpbe_getIons(Vpbe *thee, int *nion, double ionConc[MAXION],
                    double ionRadii[MAXION], double ionQ[MAXION]);

/** @brief  Object destructor
 *  @ingroup Vpbe
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void    Vpbe_dtor(Vpbe **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vpbe
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void    Vpbe_dtor2(Vpbe *thee);

/** @brief  Calculate coulombic energy of set of charges
 *         
 *           Perform an inefficient double sum to calculate the Coulombic
 *           energy of a set of charges in a homogeneous dielectric (with
 *           permittivity equal to the protein interior) and zero ionic
 *           strength.  Result is returned in units of k_B T.  The sum can be
 *           restriction to charges present in simplices of specified color
 *           (pcolor); if (color == -1) no restrictions are used.
 *
 *  @ingroup Vpbe
 *  @author  Nathan Baker
 *  @param   thee Vpbe object
 *  @return  Coulombic energy in units of \f$k_B T\f$.
 */
VEXTERNC double  Vpbe_getCoulombEnergy1(Vpbe *thee);

/** @brief   Return the memory used by this structure (and its contents)
 *           in bytes
 *  @ingroup Vpbe
 *  @author  Nathan Baker
 *  @param   thee  Vpbe object
 *  @return  The memory used by this structure and its contents in bytes
 */
VEXTERNC int     Vpbe_memChk(Vpbe *thee);

#endif /* ifndef _VPBE_H_ */
