/** @defgroup Vpmg Vpmg class
 *  @brief  A wrapper for Mike Holst's PMG multigrid code.  
 *  @note   Many of the routines and macros are borrowed from the main.c driver
 *          (written by Mike Holst) provided with the PMG code.
 */

/**
 *  @file     vpmg.h
 *  @ingroup  Vpmg
 *  @brief    Contains declarations for class Vpmg
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


#ifndef _VPMG_H_
#define _VPMG_H_

/* Generic headers */
#include "maloc/maloc.h"
#include "apbs/vhal.h"

/* Headers specific to this file */
#include "apbs/vpmgp.h"
#include "apbs/vacc.h"
#include "apbs/vcap.h"
#include "apbs/vpbe.h"
#include "apbs/vgrid.h"

/** @def VPMGMAXPART The maximum number of partitions the
 *                   mesh can be divided into 
 *  @ingroup Vpmg
 */
#define VPMGMAXPART 2000  

/** 
 *  @struct  Vpmg
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vpmg class/module
 *      
 *  Many of the routines and macros are borrowed from the main.c driver 
 *  (written by Mike Holst) provided with the PMG code.
 *     
 */
struct Vpmg {

  Vmem *vmem;                    /**< Memory management object for this class */
  Vpmgp *pmgp;                   /**< Parameters */
  Vpbe *pbe;                     /**< Information about the PBE system */

  int *iparm;                    /**< Passing int parameters to FORTRAN */
  double *rparm;                 /**< Passing real parameters to FORTRAN */
  int *iwork;                    /**< Work array */
  double *rwork;                 /**< Work array */
  double *a1cf;                  /**< Operator coefficient values (a11) */
  double *a2cf;                  /**< Operator coefficient values (a22) */
  double *a3cf;                  /**< Operator coefficient values (a33) */
  double *ccf;                   /**< Helmholtz term */
  double *fcf;                   /**< Right-hand side */
  double *tcf;                   /**< True solution */
  double *u;                     /**< Solution */
  double *xf;                    /**< Mesh point x coordinates */
  double *yf;                    /**< Mesh point y coordinates */
  double *zf;                    /**< Mesh point z coordinates */
  double *gxcf;                  /**< Boundary conditions for x faces */
  double *gycf;                  /**< Boundary conditions for y faces */
  double *gzcf;                  /**< Boundary conditions for z faces */
  int *pvec;                     /**< Partition mask array */
  double extDiEnergy;            /**< Stores contributions to the dielectric 
                                  *   energy from regions outside the problem
                                  *   domain */
  double extQmEnergy;            /**< Stores contributions to the mobile ion
                                  *   energy from regions outside the problem
                                  *   domain */
  double extQfEnergy;            /**< Stores contributions to the fixed charge
                                  *   energy from regions outside the problem
                                  *   domain */
  double extNpEnergy;            /**< Stores contributions to the apolar
                                  *   energy from regions outside the problem
                                  *   domain */
  Vsurf_Meth surfMeth;           /**< Surface definition method */
  double splineWin;              /**< Spline window parm for surf defs */
  Vchrg_Meth chargeMeth;         /**< Charge discretization method */
  int filled;                    /**< Indicates whether Vpmg_fillco has been
                                  * called */
  int useDielXMap;               /**< Indicates whether Vpmg_fillco was called
                                  * with an external x-shifted dielectric map */
  Vgrid *dielXMap;               /**< External x-shifted dielectric map */
  int useDielYMap;               /**< Indicates whether Vpmg_fillco was called
                                  * with an external y-shifted dielectric map */
  Vgrid *dielYMap;               /**< External y-shifted dielectric map */
  int useDielZMap;               /**< Indicates whether Vpmg_fillco was called
                                  * with an external z-shifted dielectric map */
  Vgrid *dielZMap;               /**< External z-shifted dielectric map */
  int useKappaMap;               /**< Indicates whether Vpmg_fillco was called
                                  * with an external kappa map */
  Vgrid *kappaMap;               /**< External kappa map */
  int useChargeMap;              /**< Indicates whether Vpmg_fillco was called
                                  * with an external charge distribution map */
  Vgrid *chargeMap;              /**< External charge distribution map */
};

/** @typedef Vpmg
 *  @ingroup Vpmg
 *  @brief   Declaration of the Vpmg class as the Vpmg structure
 */
typedef struct Vpmg Vpmg;

/* /////////////////////////////////////////////////////////////////////////
/// Inlineable methods
//////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VPMG)

    /** @brief   Return the memory used by this structure (and its contents) 
     *           in bytes
     *  @ingroup Vpmg
     *  @author  Nathan Baker
     *  @param   thee  Vpmg object
     *  @return  The memory used by this structure and its contents in bytes
     */
    VEXTERNC int Vpmg_memChk(Vpmg *thee);

#else /* if defined(VINLINE_VPMG) */

#   define Vpmg_memChk(thee) (Vmem_bytes((thee)->vmem))

#endif /* if !defined(VINLINE_VPMG) */

/* /////////////////////////////////////////////////////////////////////////
/// Non-inlineable methods
//////////////////////////////////////////////////////////////////////////// */
/** @brief   Constructor for the Vpmg class (no focusing)
 *  @author  Nathan Baker
 *  @ingroup Vpmg
 *  @param   parms  PMG parameter object
 *  @param   pbe    Object containing PBE-specific variables
 *  @returns Pointer to newly allocated Vpmg object 
 */
VEXTERNC Vpmg* Vpmg_ctor(Vpmgp *parms, Vpbe *pbe);

/** @brief   FORTRAN stub constructor for the Vpmg class (no focusing)
 *  @author  Nathan Baker
 *  @ingroup Vpmg
 *  @param   thee   Pointer to newly-allocated Vpmg object
 *  @param   parms  PMG parameter object
 *  @param   pbe    Object containing PBE-specific variables
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_ctor2(Vpmg *thee, Vpmgp *parms, Vpbe *pbe);

/** @brief   Constructor for the Vpmg class (with focusing)
 *  @note    Construct the PMG object by focusing.  In other words, use the
 *           solution from the passed Vpmg object to set the boundary
 *           conditions for the new Vpmg object.  IN THE PROCESS, THE OLD VPMG
 *           OBJECT IS DESTROYED.  The solver parameters specifed by the passed
 *           Vpmgp object and the equation data from the Vpbe object are also
 *           used.
 *  @author  Nathan Baker
 *  @ingroup Vpmg
 *  @param   parms  PMG parameter object for new mesh
 *  @param   pbe    PBE parameter object
 *  @param   pmgOLD Old Vpmg object to use for setting boundary conditions
 *  @param   energyFlag  
 *             \li 0:  Don't calculate any energy contribution from
 *                   outside focusing area 
 *             \li 1:  Calculate total energy contribution from outside 
 *                   focusing area
 *             \li 2:  Calculate energy component contributions
 *  @returns Pointer to the newly allocated Vpmg object
 */
VEXTERNC Vpmg* Vpmg_ctorFocus(Vpmgp *parms, Vpbe *pbe, Vpmg *pmgOLD,
  int energyFlag);

/** @brief   FORTRAN stub bonstructor for the Vpmg class (with focusing)
 *  @note    Construct the PMG object by focusing.  In other words, use the
 *           solution from the passed Vpmg object to set the boundary
 *           conditions for the new Vpmg object.  IN THE PROCESS, THE OLD VPMG
 *           OBJECT IS DESTROYED.  The solver parameters specifed by the passed
 *           Vpmgp object and the equation data from the Vpbe object are also
 *           used.
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @param   thee   Pointer to newly allocated Vpmg object
 *  @param   parms  PMG parameter object for new mesh
 *  @param   pbe    PBE parameter object 
 *  @param   pmgOLD Old Vpmg object to use for setting boundary conditions
 *  @param   energyFlag  
 *             \li 0:  Don't calculate any energy contribution from
 *                   outside focusing area 
 *             \li 1:  Calculate total energy contribution from outside
 *                   focusing area
 *             \li 2:  Calculate energy component contributions
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_ctor2Focus(Vpmg *thee, Vpmgp *parms, Vpbe *pbe, Vpmg *pmgOLD,
  int energyFlag);

/** @brief   Object destructor
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void Vpmg_dtor(Vpmg **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void Vpmg_dtor2(Vpmg *thee);

/** @brief   Fill the coefficient arrays prior to solving the equation
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @param   thee  Vpmg object
 *  @param   surfMeth  Surface discretization method
 *  @param   splineWin    Spline window (for surfMeth = VSM_SPLINE)
 *  @param   chargeMeth   Charge discretization method
 *  @param   useDielXMap  Specifies whether to use (1) or ignore (0) the 
 *                        dielXMap arguement
 *  @param   dielXMap     Pointer to a Vgrid object containing an external
 *                        x-shifted dielectric map.  Can be VNULL if 
 *                        useDielMap is 0.
 *  @param   useDielYMap  Specifies whether to use (1) or ignore (0) the 
 *                        dielYMap arguement
 *  @param   dielYMap     Pointer to a Vgrid object containing an external
 *                        y-shifted dielectric map.  Can be VNULL if 
 *                        useDielMap is 0.
 *  @param   useDielZMap  Specifies whether to use (1) or ignore (0) the 
 *                        dielZMap arguement
 *  @param   dielZMap     Pointer to a Vgrid object containing an external
 *                        z-shifted dielectric map.  Can be VNULL if 
 *                        useDielMap is 0.
 *  @param   useKappaMap  Specifies whether to use (1) or ignore (0) the 
 *                        kappaMap arguement
 *  @param   kappaMap     Pointer to a Vgrid object containing an external
 *                        kappa map.  Can be VNULL if useKappaMap is 0.
 *  @param   useChargeMap Specifies whether to use (1) or ignore (0) the 
 *                        chargeMap arguement
 *  @param   chargeMap    Pointer to a Vgrid object containing an external
 *                        charge distribution map.  Can be VNULL if
 *                        useChargeMap is 0.
 */
VEXTERNC void Vpmg_fillco(Vpmg *thee, 
  Vsurf_Meth surfMeth,      double splineWin,  Vchrg_Meth chargeMeth,
  int useDielXMap,   Vgrid *dielXMap, 
  int useDielYMap,   Vgrid *dielYMap, 
  int useDielZMap,   Vgrid *dielZMap, 
  int useKappaMap,   Vgrid *kappaMap,
  int useChargeMap,  Vgrid *chargeMap);

/** @brief   Solve the PBE using PMG
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @param   thee  Vpmg object
 */
VEXTERNC void Vpmg_solve(Vpmg *thee);

/** @brief   Get the total electrostatic energy.
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @param   thee   Vpmg object
 *  @param   extFlag If this was a focused calculation, then it is possible
 *                   to include the energy contributions from the outside
 *                   the focused domain.  This should be on (=1) for
 *                   sequential focusing calculations and off (=0) for
 *                   parallel calculations.
 *  @returns The electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double Vpmg_energy(Vpmg *thee, int extFlag);

/** @brief   Get the "fixed charge" contribution to the electrostatic energy 
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the fixed charges
 *           with the potential: \f[ G = \sum_i q_i u(r_i) \f]
 *           and return the result in units of \f$k_B T\f$.  Clearly, no
 *           self-interaction terms are removed.  A factor a 1/2 has to be
 *           included to convert this to a real energy.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @param   thee   Vpmg object
 *  @param   extFlag If this was a focused calculation, then it is possible
 *                   to include the energy contributions from the outside
 *                   the focused domain.  This should be on (=1) for
 *                   sequential focusing calculations and off (=0) for
 *                   parallel calculations.
 *  @returns The fixed charge electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double Vpmg_qfEnergy(Vpmg *thee, int extFlag);

/** @brief   Get the per-atom "fixed charge" contribution to the electrostatic
 *           energy
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the fixed charges
 *           with the potential: \f[ G = q u(r), \f] where \f$q\f$ is the
 *           charge and \f$r\f$ is the location of the atom of interest.  The
 *           result is returned in units of \f$k_B T\f$.  Clearly, no
 *           self-interaction terms are removed.  A factor a 1/2 has to be
 *           included to convert this to a real energy.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @param   thee   Vpmg object
 *  @param   atom   The atom for which to calculate the energy.
 *  @returns The fixed charge electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double Vpmg_qfAtomEnergy(Vpmg *thee, Vatom *atom);

/** @brief Get the "mobile charge" contribution to the electrostatic energy.
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the mobile charges 
 *           with the potential: 
 *              \f[ G = \frac{1}{4 I_s} \sum_i c_i q_i^2 \int
 *              \overline{\kappa}^2(x) e^{-q_i u(x)} dx \f]
 *           for the NPBE and
 *              \f[ G = \frac{1}{2} \int \overline{\kappa}^2(x) u^2(x) dx \f]
 *           for the LPBE.  Here \f$i\f$ denotes the counterion species, 
 *           \f$I_s\f$ is the bulk ionic strength, \f$\overline{\kappa}^2(x)\f$
 *           is the modified Debye-Huckel parameter, \f$c_i\f$ is the 
 *           concentration of species \f$i\f$, \f$q_i\f$ is the charge of
 *           species \f$i\f$, and \f$u(x)\f$ is the dimensionless electrostatic
 *           potential.  The energy is scaled to units of \f$k_b T\f$.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @param   thee   Vpmg object
 *  @param   extFlag If this was a focused calculation, then it is possible
 *                   to include the energy contributions from the outside
 *                   the focused domain.  This should be on (=1) for
 *                   sequential focusing calculations and off (=0) for
 *                   parallel calculations.
 *  @returns The mobile charge electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double Vpmg_qmEnergy(Vpmg *thee, int extFlag);

/** @brief Get the "polarization" contribution to the electrostatic energy.
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the mobile charges 
 *           with the potential: 
 *              \f[ G = \frac{1}{2} \int \epsilon (\nabla u)^2 dx \f]
 *           where \f$\epsilon\f$ is the dielectric parameter and \f$u(x)\f$ is
 *           the dimensionless electrostatic potential.  The energy is scaled
 *           to units of \f$k_b T\f$.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @param   thee   Vpmg object
 *  @param   extFlag If this was a focused calculation, then it is possible
 *                   to include the energy contributions from the outside
 *                   the focused domain.  This should be on (=1) for
 *                   sequential focusing calculations and off (=0) for
 *                   parallel calculations.
 *  @returns The polarization electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double Vpmg_dielEnergy(Vpmg *thee, int extFlag);

/** @brief Get the "apolar" energy
 *
 *           Using the dielectric map at the finest mesh level, calculate the
 *           surface area in a manner consistent with the force evaluation
 *           routines of Im et al (see Vpmg_dbnpForce and Vpmg_dielGradNorm):
 *              \f[ A = \frac{1}{\epsilon_s-\epsilon_p} \int \| \nabla
 *              \epsilon \| dx \f]
 *           where \f$\epsilon\f$ is the dielectric parameter, \f$\epsilon_s\f$
 *           is the dielectric constant for the solvent and \f$\epsilon_p\f$ is
 *           the dielectric constant for the protein.  The apolar energy is
 *           then,
 *              \f[G_{np} = \gamma S \f]
 *           where \f$\gamma\f$ is the apolar coefficient set in Vpbe (see
 *           Vpbe_ctor).  The energy is returned in units of \f$k_b T\f$.
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    I personally feel that this routine should not find its way into
 *           the main APBS driver.  In this case, the apolar energy is
 *           calculated in a manner consistent with the force evaluation, but
 *           it is not the only possible apolar energy definition...
 *           The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @param   thee   Vpmg object
 *  @param   extFlag If this was a focused calculation, then it is possible
 *                   to include the energy contributions from the outside
 *                   the focused domain.  This should be on (=1) for
 *                   sequential focusing calculations and off (=0) for
 *                   parallel calculations.
 *  @returns The apolar energy in units of \f$k_B T\f$.
 */
VEXTERNC double Vpmg_npEnergy(Vpmg *thee, int extFlag);

/** @brief Get the integral of the gradient of the dielectric function
 *
 *           Using the dielectric map at the finest mesh level, calculate the
 *           integral of the norm of the dielectric function gradient
 *           routines of Im et al (see Vpmg_dbnpForce for reference):
 *              \f[ \int \| \nabla \epsilon \| dx \f]
 *           where \f$\epsilon\f$ is the dielectric parameter.
 *           The integral is returned in units of A^2.
 * 
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @param   thee   Vpmg object
 *  @returns The integral in units of A^2.
 */
VEXTERNC double Vpmg_dielGradNorm(Vpmg *thee);

/** @brief    Calculate the total force on the specified atom in units of
 *            \f$k_B T/\AA\f$
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    \li Using the force evaluation methods of Im et al (Roux group),
 *             Comput Phys Commun, 111, 59--75 (1998).  However, this gives the
 *             whole (self-interactions included) force -- reaction field
 *             forces will have to be calculated at higher level.
 *           \li No contributions are made from higher levels of focusing.
 *           \li This is currently implemented in a very inefficient fashion
 *             becuase I'm not sure which of the PMG coefficient arrays can be
 *           re-used and which are overwritten by PMG.
 *  @param   thee  Vpmg object
 *  @param   force 3*double space to hold the force in units of \f$k_B T/\AA\f$
 *  @param   atomID  Valist ID of desired atom
 *  @param   srfm    Surface definition method
 *  @param   chgm    Charge discretization method
 */
VEXTERNC void Vpmg_force(Vpmg *thee, double *force, int atomID, 
  Vsurf_Meth srfm, Vchrg_Meth chgm);

/** @brief    Calculate the "charge-field" force on the specified atom in units
 *           of \f$k_B T/\AA\f$
 * @ingroup  Vpmg
 * @author   Nathan Baker
 * @note     \li Using the force evaluation methods of Im et al (Roux group),
 *             Comput Phys Commun, 111, 59--75 (1998).  However, this gives the
 *             whole (self-interactions included) force -- reaction field 
 *             forces will have to be calculated at higher level.
 *           \li No contributions are made from higher levels of focusing. 
 *           \li This is currently implemented in a very inefficient fashion
 *             becuase I'm not sure which of the PMG coefficient arrays can be
 *           re-used and which are overwritten by PMG.
 * @param    thee  Vpmg object
 * @param    force 3*double space to hold the force in units of \f$k_B T/\AA\f$
 * @param    atomID  Valist ID of desired atom
 * @param    chgm    Charge discretization method
 */
VEXTERNC void Vpmg_qfForce(Vpmg *thee, double *force, int atomID, 
  Vchrg_Meth chgm);

/** @brief   Calculate the dielectric boundary and apolar forces on the
 *           specified atom in units of \f$k_B T/\AA\f$
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    \li Using the force evaluation methods of Im et al (Roux group),
 *             Comput Phys Commun, 111, 59--75 (1998).  However, this gives the
 *             whole (self-interactions included) force -- reaction field 
 *             forces will have to be calculated at higher level.
 *           \li No contributions are made from higher levels of focusing. 
 *           \li This is currently implemented in a very inefficient fashion
 *             becuase I'm not sure which of the PMG coefficient arrays can be
 *           re-used and which are overwritten by PMG.
 *  @param   thee  Vpmg object
 *  @param   dbForce 3*double space to hold the dielectric boudnary force in
 *           units of \f$k_B T/\AA\f$ 
 *  @param   npForce 3*double space to hold the apolar boudnary force in
 *           units of \f$k_B T/\AA\f$ 
 *  @param   atomID  Valist ID of desired atom
 *  @param   srfm    Surface definition method
 */
VEXTERNC void Vpmg_dbnpForce(Vpmg *thee, double *dbForce, double *npForce,
  int atomID, Vsurf_Meth srfm);

/** @brief   Calculate the osmotic pressure on the specified atom in units of
 *           \f$k_B T/\AA\f$
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    \li Using the force evaluation methods of Im et al (Roux group),
 *             Comput Phys Commun, 111, 59--75 (1998).  However, this gives the
 *             whole (self-interactions included) force -- reaction field 
 *             forces will have to be calculated at higher level.
 *           \li No contributions are made from higher levels of focusing. 
 *           \li This is currently implemented in a very inefficient fashion
 *             becuase I'm not sure which of the PMG coefficient arrays can be
 *           re-used and which are overwritten by PMG.
 *  @param   thee    Vpmg object
 *  @param   force   3*double space to hold the force in units of \f$k_B
 *                   T/\AA\f$ 
 *  @param   atomID  Valist ID of desired atom
 *  @param   srfm    Surface definition method
 */
VEXTERNC void Vpmg_ibForce(Vpmg *thee, double *force, int atomID, 
  Vsurf_Meth srfm);

/** @brief   Set partition information which restricts the calculation of
 *           observables to a (rectangular) subset of the problem domain
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @param   thee          Vpmg object
 *  @param   lowerCorner   Partition lower corner
 *  @param   upperCorner   Partition upper corner
 *  @param   bflags        Whether or not a particular processor owns a face of
 *                         it's partition.  This keeps things disjoint.  0 if
 *                         the face is not owned by this partition and 1
 *                         otherwise.
 */
VEXTERNC void Vpmg_setPart(Vpmg *thee, double lowerCorner[3],
  double upperCorner[3], int bflags[6]);

/** @brief   Remove partition restrictions
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @param   thee   Vpmg object 
 */
VEXTERNC void Vpmg_unsetPart(Vpmg *thee);

/** @brief   Fill the specified array with accessibility values 
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @param   thee   Vpmg object
 *  @param   vec    An array nx*ny*nz in length to contain the accessibility
 *                  values (where nx, ny, nz are the numbers of grid points)
 *  @param   type   Which data to write
 *  @param   parm   Parameter for data definition (if needed)
 */
VEXTERNC void Vpmg_fillArray(Vpmg *thee, double *vec, Vdata_Type type, 
  double parm);

/** @brief   Print out a column-compressed sparse matrix in Harwell-Boeing
 *           format.  
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @param   thee   Vpmg object
 *  @param   path   The file to which the matrix is to be written
 *  @param   title  The title of the matrix
 *  @param   mxtype The type of REAL-valued matrix, a 3-character string of the
 *                  form "R_A", where the blank can contain:
 *                  \li S:  Symmetric matrix
 *                  \li U:  Unsymmetric matrix
 *                  \li H:  Hermitian matrix
 *                  \li Z:  Skew symmetric matrix
 *                  \li R:  Rectangular
 *  @param   flag   The operator to compress:
 *                  \li 0:  Just the Poisson operator
 *                  \li 1:  The linearization of the full Poisson-Boltzmann
 *                          operator around the current solution
 */
VEXTERNC void Vpmg_printColComp(Vpmg *thee, char path[72], char title[72],
  char mxtype[3], int flag);

#endif    /* ifndef _VPMG_H_ */
