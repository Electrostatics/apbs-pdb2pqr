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
// File:     mgparm.c
//
// Purpose:  Class MGparm: methods. 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "apbs/apbs.h"
#include "apbs/vhal.h"
#include "apbs/mgparm.h"
#include "apbs/vstring.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class MGparm: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_MGPARM)

#endif /* if !defined(VINLINE_MGPARM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class MGparm: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_ctor
//
// Author: Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC MGparm* MGparm_ctor(int type) {

    /* Set up the structure */
    MGparm *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(MGparm));
    VASSERT( thee != VNULL);
    VASSERT( MGparm_ctor2(thee, type) );

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_ctor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int MGparm_ctor2(MGparm *thee, int type) {

    int i;

    if (thee == VNULL) return 0;

    for (i=0; i<3; i++) {
        thee->dime[i] = -1;
        thee->pdime[i] = 1;
    }

    thee->parsed = 0;
    thee->type = type;

    /* *** GENERIC PARAMETERS *** */
    thee->setdime = 0;

    /* *** TYPE 0 PARAMETERS *** */
    thee->nlev = VMGNLEV;
    thee->setnlev = 1;
    thee->setgrid = 0;
    thee->setglen = 0;
    thee->setgcent = 0;  

    /* *** TYPE 1 & 2 PARAMETERS *** */
    thee->setcglen = 0;
    thee->setfglen = 0;
    thee->setcgcent = 0;
    thee->setfgcent = 0;

    /* *** TYPE 2 PARAMETERS *** */
    thee->setpdime = 0;
    thee->setcom = 0;
    thee->setofrac = 0;
    for (i=0; i<6; i++) thee->partDisjOwnSide[i] = 1;

    return 1; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_dtor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void MGparm_dtor(MGparm **thee) {
    if ((*thee) != VNULL) {
        MGparm_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(MGparm), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_dtor2
//
// Purpose:  Destroy the atom object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void MGparm_dtor2(MGparm *thee) { ; }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_check
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int MGparm_check(MGparm *thee) { 

    int rc, i, tdime[3], ti, tnlev[3], nlev;

    rc = 1;

    /* Check to see if we were even filled... */
    if (!thee->parsed) {
        Vnm_print(2, "MGparm_check:  not filled!\n");
        return 0;
    }

    /* Check generic settings */
    if (!thee->setdime) {
        Vnm_print(2, "MGparm_check:  DIME not set!\n");
        rc = 0;
    }

    /* Check sequential manual settings */
    if (thee->type == 0) {
        if ((!thee->setgrid) && (!thee->setglen)) {
            Vnm_print(2, "MGparm_check:  Neither GRID nor GLEN set!\n");
            rc = 0;
        }
        if ((thee->setgrid) && (thee->setglen)) {
            Vnm_print(2, "MGparm_check:  Both GRID and GLEN set!\n");
            rc = 0;
        }
        if (!thee->setgcent) {
            Vnm_print(2, "MGparm_check:  GCENT not set!\n");
            rc = 0;
        }
    }
 
    /* Check sequential and parallel automatic focusing settings */
    if ((thee->type == 1) || (thee->type == 2)) {
        if (!thee->setcglen) {
            Vnm_print(2, "MGparm_check:  CGLEN not set!\n");
            rc = 0;
        }
        if (!thee->setfglen) {
            Vnm_print(2, "MGparm_check:  FGLEN not set!\n");
            rc = 0;
        }
        if (!thee->setcgcent) {
            Vnm_print(2, "MGparm_check:  CGCENT not set!\n");
            rc = 0;
        }
        if (!thee->setfgcent) {
            Vnm_print(2, "MGparm_check:  FGCENT not set!\n");
            rc = 0;
        }
    }

    /* Check parallel automatic focusing settings */
    if (thee->type == 2) {
        if (!thee->setpdime) {
            Vnm_print(2, "MGparm_check:  PDIME not set!\n");
            rc = 0;
        }
        if (!thee->setcom) {
            Vnm_print(2, "MGparm_check:  COM not set!\n");
            rc = 0;
        }
        if (!thee->setofrac) {
            Vnm_print(2, "MGparm_check:  OFRAC not set!\n");
            rc = 0;
        }
    }
 
    /* Perform a sanity check on nlev and dime, resetting values as necessary */
    if (rc == 1) {
	/* Calculate the actual number of grid points and nlev to satisfy the
	 * formula:  n = c * 2^(l+1) + 1, where n is the number of grid points,
	 * c is an integer, and l is the number of levels */
        for (i=0; i<3; i++) {
            tdime[i] = thee->dime[i];
            ti = tdime[i] - 1;
            tnlev[i] = 0;
            /* Find the maximum number of times this dimension can be divided by
             * two */
            while (1) {
                if (VODD(ti)) break;
                (tnlev[i])++;
                ti = (int)ceil(0.5*ti);
            }
            (tnlev[i])--;
	    /* We'd like to have at least VMGNLEV levels in the multigrid
	     * hierarchy.  This means that the dimension needs to be c*2^VMGNLEV
	     * + 1, where c is an integer. */
            if (tnlev[i] < VMGNLEV) {
                Vnm_print(2, "NOsh:  Bad dime[%d]  = %d (%d nlev)!\n",
                  i, tdime[i], tnlev[i]);
                ti = (int)(tdime[i]/VPOW(2.,(VMGNLEV+1)));
                if (ti < 1) ti = 1;
                tdime[i] = ti*VPOW(2.,(VMGNLEV+1)) + 1;
                tnlev[i] = 4;
                Vnm_print(2, "NOsh:  Reset dime[%d] to %d and (nlev = %d).\n", 
                  i, tdime[i], VMGNLEV);
            }
        }
	/* The actual number of levels we'll be using is the smallest number of
         * possible levels in any dimensions */
        nlev = VMIN2(tnlev[0], tnlev[1]);
        nlev = VMIN2(nlev, tnlev[2]);
        /* Set the number of levels and dimensions */
        Vnm_print(0, "NOsh:  nlev = %d, dime = (%d, %d, %d)\n", nlev, tdime[0],
          tdime[1], tdime[2]);
        thee->nlev = nlev;
        for (i=0; i<3; i++) thee->dime[i] = tdime[i];
    }

    return rc;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_copy
//  
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void MGparm_copy(MGparm *thee, MGparm *parm) {

    int i;

    VASSERT(thee != VNULL);
    VASSERT(parm != VNULL);

    
    thee->type = parm->type;
    thee->parsed = parm->parsed;

    /* *** GENERIC PARAMETERS *** */
    for (i=0; i<3; i++) thee->dime[i] = parm->dime[i];
    thee->setdime = parm->setdime;

    /* *** TYPE 0 PARMS *** */
    thee->nlev = parm->nlev;
    thee->setnlev = parm->setnlev;
    for (i=0; i<3; i++) thee->grid[i] = parm->grid[i];
    thee->setgrid = parm->setgrid;
    for (i=0; i<3; i++) thee->glen[i] = parm->glen[i];
    thee->setglen = parm->setglen;
    thee->cmeth = parm->cmeth;
    for (i=0; i<3; i++) thee->center[i] = parm->center[i];
    thee->setgcent = parm->setgcent;
    thee->centmol = parm->centmol;

    /* *** TYPE 1 & 2 PARMS *** */
    for (i=0; i<3; i++) thee->cglen[i] = parm->cglen[i];
    thee->setcglen = parm->setcglen;
    for (i=0; i<3; i++) thee->fglen[i] = parm->fglen[i];
    thee->setfglen = parm->setfglen;
    thee->ccmeth = parm->ccmeth;
    for (i=0; i<3; i++) thee->ccenter[i] = parm->ccenter[i];
    thee->setcgcent = parm->setcgcent;
    thee->ccentmol = parm->ccentmol;
    thee->fcmeth = parm->fcmeth;
    for (i=0; i<3; i++) thee->fcenter[i] = parm->fcenter[i];
    thee->setfgcent = parm->setfgcent;
    thee->fcentmol = parm->fcentmol;

    /* *** TYPE 2 PARMS *** */
    for (i=0; i<3; i++) 
      thee->partDisjCenterShift[i] = parm->partDisjCenterShift[i];
    for (i=0; i<3; i++) 
      thee->partDisjLength[i] = parm->partDisjLength[i];
    for (i=0; i<3; i++) 
      thee->partDisjOwnSide[i] = parm->partDisjOwnSide[i];
    for (i=0; i<3; i++) 
      thee->partOlapCenterShift[i] = parm->partOlapCenterShift[i];
    for (i=0; i<3; i++) 
      thee->partOlapLength[i] = parm->partOlapLength[i];
    for (i=0; i<3; i++) thee->pdime[i] = parm->pdime[i];
    thee->setpdime = parm->setpdime;
    thee->com = parm->com;
    thee->setcom = parm->setcom;
    thee->ofrac = parm->ofrac;
    thee->setofrac = parm->setofrac;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGparm_parseToken
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int MGparm_parseToken(MGparm *thee, char tok[VMAX_BUFSIZE], 
  Vio *sock) {

    int ti;
    double tf;

    if (thee == VNULL) {
        Vnm_print(2, "parseMG:  got NULL thee!\n"); 
        return -1;
    }
    if (sock == VNULL) {
        Vnm_print(2, "parseMG:  got NULL socket!\n"); 
        return -1;
    }

    if (Vstring_strcasecmp(tok, "dime") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0){
            Vnm_print(2, "parseMG:  Read non-integer (%s) while parsing DIME \
keyword!\n", tok);
            return -1;
        } else thee->dime[0] = ti;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing DIME \
keyword!\n", tok);
            return -1;
        } else thee->dime[1] = ti;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing DIME \
keyword!\n", tok);
            return -1;
        } else thee->dime[2] = ti;
        thee->setdime = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "nlev") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing NLEV \
keyword!\n", tok);
            return -1;
        } else thee->nlev = ti;
        thee->setnlev = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "grid") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GRID \
keyword!\n", tok);
            return -1;
        } else thee->grid[0] = tf;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GRID \
keyword!\n", tok);
            return -1;
        } else thee->grid[1] = tf;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GRID \
keyword!\n", tok);
            return -1;
        } else thee->grid[2] = tf;
        thee->setgrid = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "glen") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN \
keyword!\n", tok);
            return -1;
        } else thee->glen[0] = tf;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN \
keyword!\n", tok);
            return -1;
        } else thee->glen[1] = tf;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN \
keyword!\n", tok);
            return -1;
        } else thee->glen[2] = tf;
        thee->setglen = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "gcent") == 0) {
        /* If the next token isn't a float, it probably means we want to
         * center on a molecule */
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            if (Vstring_strcasecmp(tok, "mol") == 0) {
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%d", &ti) == 0) {
                    Vnm_print(2, "NOsh:  Read non-int (%s) while parsing \
GCENT MOL keyword!\n", tok);
                    return -1;
                } else {
                    thee->cmeth = 1;
                    thee->centmol = ti;
                }
            } else {
                Vnm_print(2, "NOsh:  Unexpected keyword (%s) while parsing \
GCENT!\n", tok);
                return -1;
            }
        } else {
            thee->center[0] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing \
GCENT keyword!\n", tok);
                return -1;
            }
            thee->center[1] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing \
GCENT keyword!\n", tok);
                return -1;
            } 
            thee->center[2] = tf;
        }
        thee->setgcent = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "cglen") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing CGLEN \
keyword!\n", tok);
            return -1;
        } else thee->cglen[0] = tf;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing CGLEN \
keyword!\n", tok);
            return -1;
        } else thee->cglen[1] = tf;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing CGLEN \
keyword!\n", tok);
            return -1;
        } else thee->cglen[2] = tf;
        thee->setcglen = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "fglen") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing FGLEN \
keyword!\n", tok);
            return -1;
        } else thee->fglen[0] = tf;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing FGLEN \
keyword!\n", tok);
            return -1;
        } else thee->fglen[1] = tf;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing FGLEN \
keyword!\n", tok);
            return -1;
        } else thee->fglen[2] = tf;
        thee->setfglen = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "cgcent") == 0) {
        /* If the next token isn't a float, it probably means we want to
         * center on a molecule */
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            if (Vstring_strcasecmp(tok, "mol") == 0) {
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%d", &ti) == 0) {
                    Vnm_print(2, "NOsh:  Read non-int (%s) while parsing \
CGCENT MOL keyword!\n", tok);
                    return -1;
                } else {
                    thee->ccmeth = 1;
                    thee->ccentmol = ti;
                }
            } else {
                Vnm_print(2, "NOsh:  Unexpected keyword (%s) while parsing \
CGCENT!\n", tok);
                return -1;
                }
        } else {
            thee->ccenter[0] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing \
CGCENT keyword!\n", tok);
                return -1;
            }
            thee->ccenter[1] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing \
CGCENT keyword!\n", tok);
                return -1;
            }
            thee->ccenter[2] = tf;
        }
        thee->setcgcent = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "fgcent") == 0) {
        /* If the next token isn't a float, it probably means we want to
         * center on a molecule */
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            if (Vstring_strcasecmp(tok, "mol") == 0) {
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%d", &ti) == 0) {
                     Vnm_print(2, "NOsh:  Read non-int (%s) while parsing \
FGCENT MOL keyword!\n", tok);
                     return -1;
                } else {
                    thee->fcmeth = 1;
                    thee->fcentmol = ti;
                }
            } else {
                Vnm_print(2, "NOsh:  Unexpected keyword (%s) while parsing \
FGCENT!\n", tok);
                return -1;
            }
        } else {
            thee->fcenter[0] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing \
FGCENT keyword!\n", tok);
                return -1;
            }
            thee->fcenter[1] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing \
FGCENT keyword!\n", tok);
                return -1;
            }
            thee->fcenter[2] = tf;
        }
        thee->setfgcent = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "pdime") == 0) {
        /* Read the number of grid points */
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing PDIME \
keyword!\n", tok);
            return -1;
        } else {
            thee->pdime[0] = ti;
        }
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing PDIME \
keyword!\n", tok);
            return -1;
        } else {
            thee->pdime[1] = ti;
        }
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing PDIME \
keyword!\n", tok);
            return -1;
        } else {
            thee->pdime[2] = ti;
        }
        thee->setpdime = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "ofrac") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing OFRAC \
keyword!\n", tok);
            return -1;
        }
        thee->ofrac = tf;
        thee->setofrac = 1;
        return 1;
    }


    return 0;

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;

}
