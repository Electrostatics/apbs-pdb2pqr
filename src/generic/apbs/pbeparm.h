/* ///////////////////////////////////////////////////////////////////////////
/// APBS -- Adaptive Poisson-Boltzmann Solver
///
///  Nathan A. Baker (nbaker@wasabi.ucsd.edu)
///  Dept. of Chemistry and Biochemistry
///  Dept. of Mathematics, Scientific Computing Group
///  University of California, San Diego 
///
///  Additional contributing authors listed in the code documentation.
///
/// Copyright © 1999. The Regents of the University of California (Regents).
/// All Rights Reserved. 
/// 
/// Permission to use, copy, modify, and distribute this software and its
/// documentation for educational, research, and not-for-profit purposes,
/// without fee and without a signed licensing agreement, is hereby granted,
/// provided that the above copyright notice, this paragraph and the
/// following two paragraphs appear in all copies, modifications, and
/// distributions.
/// 
/// IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
/// SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
/// ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
/// REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
/// 
/// REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
/// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
/// PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
/// ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
/// TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
/// MODIFICATIONS. 
//////////////////////////////////////////////////////////////////////////// 
/// rcsid="$Id$"
//////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// File:     pbeparse.h    
//
// Purpose:  A set of useful parameters and parsing methods for a generic PBE
//           calculation
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _PBEPARM_H_
#define _PBEPARM_H_

#include "apbs/apbs.h"
#include "maloc/maloc.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class PBEparm: Definition
//
// If you add/change something here, it must be added/changed in PBEparm_copy
/////////////////////////////////////////////////////////////////////////// */
typedef struct PBEparm {

    int molid;                 /* Molecule ID to perform calculation on */
    int setmolid;
    int nonlin;                /* 0 => LPBE, 1 => NPBE */
    int setnonlin;
    int bcfl;                  /* Boundary condition: 0 => zero, 1 => single
                                * Debye-Huckel sphere, 2 => multiple Debye-
                                * Huckel spheres, 4 => focusing */
    int setbcfl;
    int nion;                  /* Number of counterion species */
    int setnion;
    double ionq[MAXION];  /* Counterion charges (in e) */
    double ionc[MAXION];  /* Counterion concentrations (in M) */
    double ionr[MAXION];  /* Counterion radii (in A) */
    int setion[MAXION];
    double pdie;               /* Solute dielectric */
    int setpdie;
    double sdie;               /* Solvent dielectric */
    int setsdie;
    int srfm;                  /* Surface calculation method
				* 0 => Mol surface for epsilon; inflated VdW
				* for kappa; no smoothing 
                                * 1 => As 0 with harmoic average
                                * smoothing
                                * 2 => Cubic spline */
    int setsrfm;
    double srad;               /* Solvent radius */
    int setsrad;
    double swin;               /* Cubic spline window */
    int setswin; 
    double temp;               /* Temperarture (in K) */
    int settemp;
    double gamma;              /* Surface tension for apolar energies/forces
                                * (in kJ/mol/A^2) */
    int setgamma;
    int calcenergy;            /* Energy calculation
				* 0 => don't calculate out energy
                                * 1 => calculate total energy 
                                * 2 => calculate total energy and all energy
                                *      components*/
    int setcalcenergy;      
    int calcforce;             /* Atomic forces I/O 
                                * 0 => don't calculate forces
                                * 1 => calculate net forces on molecule
                                * 2 => calculate atom-level forces */
    int setcalcforce;       
    int writepot;              /* 0 => no, 1 => yes */
    int setwritepot; 
    char writepotstem[VMAX_ARGLEN];
                               /* File stem to write pot */
    int writepotfmt;           /* Potential file formats: 0 => dx, 1 => avs, 
                                * 2 => UHBD */
    int writeacc;              /* 0 => no, 1 => yes */
    int setwriteacc; 
    char writeaccstem[VMAX_ARGLEN];    
                               /* File stem to write pot */
    int writeaccfmt;           /* Potential file formats: 0 => dx, 1 => avs, 
                                * 2 => UHBD */

    int parsed;                /* Has this been filled with anything other than
                                * the default values? */
} PBEparm;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (mcsh.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC PBEparm* PBEparm_ctor();
VEXTERNC int      PBEparm_ctor2(PBEparm *thee);
VEXTERNC void     PBEparm_dtor(PBEparm **thee);
VEXTERNC void     PBEparm_dtor2(PBEparm *thee);
VEXTERNC int      PBEparm_check(PBEparm *thee);
VEXTERNC void     PBEparm_copy(PBEparm *thee, PBEparm *parm);
VEXTERNC int      PBEparm_parseToken(PBEparm *thee, char tok[VMAX_BUFSIZE],
                    Vio *sock);


#endif 

