/** @defgroup NOsh NOsh class
 *  @brief    Class for parsing for fixed format input files
 */

/**
 *  @file     nosh.h
 *  @ingroup  NOsh
 *  @brief    Contains declarations for class NOsh
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

#ifndef _NOSH_H_
#define _NOSH_H_

/** @brief Maximum number of molecules in a run 
 *  @ingroup NOsh */
#define NOSH_MAXMOL 20

/** @brief Maximum number of calculations in a run 
 *  @ingroup NOsh */
#define NOSH_MAXCALC 20

/** @brief Maximum number of PRINT statements in a run 
 *  @ingroup NOsh */
#define NOSH_MAXPRINT 20

/** @brief Maximum number of operations in a PRINT statement
 *  @ingroup NOsh */
#define NOSH_MAXPOP 20

#include "apbs/apbs.h"
#include "maloc/maloc.h"
#include "apbs/femparm.h"
#include "apbs/mgparm.h"
#include "apbs/pbeparm.h"

/**
 *  @struct  NOsh_calc
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @brief   Calculation class for use when parsing fixed format input files
 */
struct NOsh_calc {

    MGparm *mgparm;         /**< Multigrid parameters */
    FEMparm *femparm;       /**< Finite element parameters */
    PBEparm *pbeparm;       /**< Generic PBE parameters */
    int calctype;           /**< 0: multigrid, 1: FEM */

};

/** @typedef NOsh_calc
 *  @ingroup NOsh
 *  @brief   Declaration of the NOsh_calc class as the NOsh_calc structure
 */
typedef struct NOsh_calc NOsh_calc;

/**
 *  @struct  NOsh
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @brief   Class for parsing fixed format input files
 */
struct NOsh {

    NOsh_calc calc[NOSH_MAXCALC];        /**< The array of calculation objects
                                          */

    int ncalc;                           /**< The number of calculations in the
					  * calc array */
    int nelec;                           /**< The number of elec statements in
                                          * the input file */
    int ispara;                          /**< 1 => is a parallel calculation, 
                                          *  0 => is not */
    int proc_rank;                       /**< Processor rank in parallel
                                          * calculation */
    int proc_size;                       /**< Number of processors in parallel
                                          * calculation */
    int bogus;                           /**< A flag which tells routines using
                                          * NOsh that this particular NOsh is
                                          * broken -- useful for parallel
                                          * focusing calculations where the
                                          * user gave us too many processors 
                                          * (1 => ignore this NOsh; 
                                          *  0 => this NOsh is OK) */
    int elec2calc[NOSH_MAXCALC];         /**< A mapping between ELEC statements
					  * which appear in the input file and
					  * calc objects stored above.  Since
					  * we allow both normal and focused
					  * multigrid, there isn't a 1-to-1
					  * correspondence between ELEC
					  * statements and actual calcualtions.
					  * This can really confuse operations
					  * which work on specific calculations
					  * further down the road (like PRINT).
					  * Therefore this array is the initial
					  * point of entry for any
					  * calculation-specific operation.  It
					  * points to a specific entry in the
					  * calc array. */
    int nmol;                            /**< Number of molecules */
    char molpath[NOSH_MAXMOL][VMAX_ARGLEN];   /**< Paths to mol files */
    int molfmt[NOSH_MAXMOL];            /**< Mol files formats (0=>PQR) */
    int ndiel;                           /**< Number of dielectric maps */
    char dielXpath[NOSH_MAXMOL][VMAX_ARGLEN]; /**< Paths to x-shifted 
                                               * dielectric map files */
    char dielYpath[NOSH_MAXMOL][VMAX_ARGLEN]; /**< Paths to y-shifted 
                                               * dielectric map files */
    char dielZpath[NOSH_MAXMOL][VMAX_ARGLEN]; /**< Paths to z-shifted 
                                               * dielectric map files */
    int dielfmt[NOSH_MAXMOL];           /**< Dielectric maps file 
                                          * formats (0=>OpenDX) */
    int nkappa;                          /**< Number of kappa maps */
    char kappapath[NOSH_MAXMOL][VMAX_ARGLEN];   /**< Paths to kappa map
                                          * files */
    int kappafmt[NOSH_MAXMOL];          /**< Kappa maps file 
                                          * formats (0=>OpenDX) */
    int ncharge;                         /**< Number of charge maps */
    char chargepath[NOSH_MAXMOL][VMAX_ARGLEN];   /**< Paths to charge map
                                          * files */
    int chargefmt[NOSH_MAXMOL];         /**< Charge maps file 
                                          * formats (0=>OpenDX) */
    int nprint;                          /**< How many print sections? */
    int printwhat[NOSH_MAXPRINT];        /**< What do we print:
                                          * \li 0 = energy
                                          * \li 1 = force */
    int printnarg[NOSH_MAXPRINT];        /**< How many arguments in energy 
                                          * list */
    int printcalc[NOSH_MAXPRINT][NOSH_MAXPOP]; /**< ELEC id (see elec2calc) */
    int printop[NOSH_MAXPRINT][NOSH_MAXPOP];  /**< Operation id (0 = add, 1 = 
                                          * subtract) */
  int parsed;                            /**< Have we parsed an input file
                                          * yet? */

};

/** @typedef NOsh
 *  @ingroup NOsh
 *  @brief   Declaration of the NOsh class as the NOsh structure
 */ 
typedef struct NOsh NOsh;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Inlineable methods (mcsh.c)
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_NOSH)
    /** @brief    Returns path to specified molecule
     *  @ingroup  NOsh
     *  @author   Nathan Baker
     *  @param    thee Pointer to NOsh object
     *  @param    imol Molecule ID of interest
     *  @returns  Path string
     */
    VEXTERNC char* NOsh_getMolpath(NOsh *thee, int imol);

    /** @brief    Returns path to specified x-shifted dielectric map
     *  @ingroup  NOsh
     *  @author   Nathan Baker
     *  @param    thee Pointer to NOsh object
     *  @param    imap Map ID of interest
     *  @returns  Path string
     */
    VEXTERNC char* NOsh_getDielXpath(NOsh *thee, int imap);

    /** @brief    Returns path to specified y-shifted dielectric map
     *  @ingroup  NOsh
     *  @author   Nathan Baker
     *  @param    thee Pointer to NOsh object
     *  @param    imap Map ID of interest
     *  @returns  Path string
     */
    VEXTERNC char* NOsh_getDielYpath(NOsh *thee, int imap);

    /** @brief    Returns path to specified z-shifted dielectric map
     *  @ingroup  NOsh
     *  @author   Nathan Baker
     *  @param    thee Pointer to NOsh object
     *  @param    imap Map ID of interest
     *  @returns  Path string
     */
    VEXTERNC char* NOsh_getDielZpath(NOsh *thee, int imap);

    /** @brief    Returns path to specified kappa map
     *  @ingroup  NOsh
     *  @author   Nathan Baker
     *  @param    thee Pointer to NOsh object
     *  @param    imap Map ID of interest
     *  @returns  Path string
     */
    VEXTERNC char* NOsh_getKappapath(NOsh *thee, int imap);

    /** @brief    Returns path to specified charge distribution map
     *  @ingroup  NOsh
     *  @author   Nathan Baker
     *  @param    thee Pointer to NOsh object
     *  @param    imap Map ID of interest
     *  @returns  Path string
     */
    VEXTERNC char* NOsh_getChargepath(NOsh *thee, int imap);

    /** @brief    Returns specified calculation object
     *  @ingroup  NOsh
     *  @author   Nathan Baker
     *  @param    thee Pointer to NOsh object
     *  @param    icalc Calculation ID of interest
     *  @returns  Pointer to specified calculation object
     */
    VEXTERNC NOsh_calc* NOsh_getCalc(NOsh *thee, int icalc);

    /** @brief    Returns format of specified dielectric map
     *  @ingroup  NOsh
     *  @author   Nathan Baker
     *  @param    thee Pointer to NOsh object
     *  @param    imap Calculation ID of interest
     *  @returns  Format of dielectric map
     */
    VEXTERNC int NOsh_getDielfmt(NOsh *thee, int imap);

    /** @brief    Returns format of specified kappa map
     *  @ingroup  NOsh
     *  @author   Nathan Baker
     *  @param    thee Pointer to NOsh object
     *  @param    imap Calculation ID of interest
     *  @returns  Format of kappa map
     */
    VEXTERNC int NOsh_getKappafmt(NOsh *thee, int imap);

    /** @brief    Returns format of specified charge map
     *  @ingroup  NOsh
     *  @author   Nathan Baker
     *  @param    thee Pointer to NOsh object
     *  @param    imap Calculation ID of interest
     *  @returns  Format of charge map
     */
    VEXTERNC int NOsh_getChargefmt(NOsh *thee, int imap);

#else

#   define NOsh_getMolpath(thee, imol) ((thee)->molpath[(imol)])
#   define NOsh_getDielXpath(thee, imol) ((thee)->dielXpath[(imol)])
#   define NOsh_getDielYpath(thee, imol) ((thee)->dielYpath[(imol)])
#   define NOsh_getDielZpath(thee, imol) ((thee)->dielZpath[(imol)])
#   define NOsh_getKappapath(thee, imol) ((thee)->kappapath[(imol)])
#   define NOsh_getChargepath(thee, imol) ((thee)->chargepath[(imol)])
#   define NOsh_getCalc(thee, icalc) ((thee)->calc[(icalc)])
#   define NOsh_getDielfmt(thee, imap) ((thee)->dielfmt[(imap)])
#   define NOsh_getKappafmt(thee, imap) ((thee)->kappafmt[(imap)])
#   define NOsh_getChargefmt(thee, imap) ((thee)->chargefmt[(imap)])

#endif


/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (mcsh.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Return an integer ID of the observable to print (@see printwhat)
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee NOsh object to use
 *  @param   iprint ID of PRINT statement
 *  @returns An integer ID of the observable to print (@see printwhat)
 */
VEXTERNC int NOsh_printWhat(NOsh *thee, int iprint);

/** @brief   Return an integer mapping of an ELEC statement to a calculation ID
 *           (@see elec2calc)
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee NOsh object to use
 *  @param   icalc ID of CALC statement
 *  @returns An integer mapping of an ELEC statement to a calculation ID
 *           (@see elec2calc)
 */
VEXTERNC int NOsh_elec2calc(NOsh *thee, int iprint);

/** @brief   Return number of arguments to PRINT statement (@see printnarg)
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee NOsh object to use
 *  @param   iprint ID of PRINT statement
 *  @returns Number of arguments to PRINT statement (@see printnarg)
 */
VEXTERNC int NOsh_printNarg(NOsh *thee, int iprint);

/** @brief   Return integer ID for specified operation (@see printop)
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee NOsh object to use
 *  @param   iprint ID of PRINT statement
 *  @param   iarg ID of operation in PRINT statement
 *  @returns Integer ID for specified operation (@see printop)
 */ 
VEXTERNC int NOsh_printOp(NOsh *thee, int iprint, int iarg);

/** @brief   Return calculation ID for specified PRINT statement 
 *           (@see printcalc)
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee NOsh object to use
 *  @param   iprint ID of PRINT statement
 *  @param   iarg ID of operation in PRINT statement
 *  @returns Calculation ID for specified PRINT statement 
 *           (@see printcalc)
 */
VEXTERNC int NOsh_printCalc(NOsh *thee, int iprint, int iarg);

/** @brief   Construct NOsh
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   rank   Rank of current processor in parallel calculation (0 if not
 *                  parallel)
 *  @param   size   Number of processors in parallel calculation (1 if not
 *                  parallel)
 *  @returns Newly allocated and initialized NOsh object
 */
VEXTERNC NOsh* NOsh_ctor(int rank, int size);

/** @brief   FORTRAN stub to construct NOsh
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee  Space for NOsh objet
 *  @param   rank   Rank of current processor in parallel calculation (0 if not
 *                  parallel)
 *  @param   size   Number of processors in parallel calculation (1 if not
 *                  parallel)
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int   NOsh_ctor2(NOsh *thee, int rank, int size);

/** @brief   Object destructor
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location of NOsh object
 */
VEXTERNC void  NOsh_dtor(NOsh **thee);

/** @brief   FORTRAN stub for object destructor
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee  Pointer to NOsh object
 */
VEXTERNC void  NOsh_dtor2(NOsh *thee);

/** @brief   Parse an input file from a socket
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee  Pointer to NOsh object
 *  @param   sock  Stream of tokens to parse
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int   NOsh_parse(NOsh *thee, Vio *sock);

/** @brief   Parse an input file only from a file
 *  @note    Included for SWIG wrapper compatibility
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee      Pointer to NOsh object
 *  @param   filename  Name/path of readable file
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int   NOsh_parseFile(NOsh *thee, char *filename);

/** @brief   Setup NOsh, MGparm, and PBEparm objects for a MG-MANUAL ELEC
 *           calculation
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee    Pointer to NOsh object
 *  @param   mgparm  Pointer to MGparm object containing basic parameters for
 *                   run; these will be supplemented and modified (for
 *                   consistency, etc.) by this function
 *  @param   pbeparm Pointer to PBEparm object containing basic parameters for
 *                   run; these will be supplemented and modified (for
 *                   consistency, etc.) by this function
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int NOsh_setupMGMANUAL(NOsh *thee, MGparm *mgparm, PBEparm *pbeparm);

/** @brief   Setup NOsh, MGparm, and PBEparm objects for a MG-AUTO ELEC
 *           calculation
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee    Pointer to NOsh object
 *  @param   mgparm  Pointer to base MGparm object containing parameters for
 *                   the run; these will be used to construct the MGparm
 *                   objects inside the NOsh object (for actual use by the code
 *                   during calculations)
 *  @param   pbeparm Pointer to PBEparm object containing basic parameters for
 *                   run; these will be supplemented and modified (for
 *                   consistency, etc.) by this function
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int NOsh_setupMGAUTO(NOsh *thee, MGparm *mgparm, PBEparm *pbeparm);

/** @brief   Setup NOsh, MGparm, and PBEparm objects for a MG-PARA ELEC
 *           calculation
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @param   thee    Pointer to NOsh object
 *  @param   mgparm  Pointer to base MGparm object containing parameters for
 *                   the run; these will be used to construct the MGparm
 *                   objects inside the NOsh object (for actual use by the code
 *                   during calculations)
 *  @param   pbeparm Pointer to PBEparm object containing basic parameters for
 *                   run; these will be supplemented and modified (for
 *                   consistency, etc.) by this function
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int NOsh_setupMGPARA(NOsh *thee, MGparm *mgparm, PBEparm *pbeparm);


#endif 

