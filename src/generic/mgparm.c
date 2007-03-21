/**
 *  @file    mgparm.c
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @brief   Class MGparm methods
 *  @version $Id$
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
 * Copyright (c) 2002-2007.  Washington University in St. Louis.
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

#include "apbscfg.h"
#include "apbs/apbs.h"
#include "apbs/vhal.h"
#include "apbs/mgparm.h"
#include "apbs/vstring.h"

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_MGPARM)

#endif /* if !defined(VINLINE_MGPARM) */

VPUBLIC void MGparm_setCenterX(MGparm *thee, double x) {
    VASSERT(thee != VNULL);
    thee->center[0] = x;
}
VPUBLIC void MGparm_setCenterY(MGparm *thee, double y) {
    VASSERT(thee != VNULL);
    thee->center[1] = y;
}
VPUBLIC void MGparm_setCenterZ(MGparm *thee, double z) {
    VASSERT(thee != VNULL);
    thee->center[2] = z;
}
VPUBLIC double MGparm_getCenterX(MGparm *thee) {
    VASSERT(thee != VNULL);
    return thee->center[0];
}
VPUBLIC double MGparm_getCenterY(MGparm *thee) {
    VASSERT(thee != VNULL);
    return thee->center[1];
}
VPUBLIC double MGparm_getCenterZ(MGparm *thee) {
    VASSERT(thee != VNULL);
    return thee->center[2];
}
VPUBLIC int MGparm_getNx(MGparm *thee) {
    VASSERT(thee != VNULL);
    return thee->dime[0];
}
VPUBLIC int MGparm_getNy(MGparm *thee) {
    VASSERT(thee != VNULL);
    return thee->dime[1];
}
VPUBLIC int MGparm_getNz(MGparm *thee) {
    VASSERT(thee != VNULL);
    return thee->dime[2];
}
VPUBLIC double MGparm_getHx(MGparm *thee) {
    VASSERT(thee != VNULL);
    return thee->grid[0];
}
VPUBLIC double MGparm_getHy(MGparm *thee) {
    VASSERT(thee != VNULL);
    return thee->grid[1];
}
VPUBLIC double MGparm_getHz(MGparm *thee) {
    VASSERT(thee != VNULL);
    return thee->grid[2];
}

VPUBLIC MGparm* MGparm_ctor(MGparm_CalcType type) {

    /* Set up the structure */
    MGparm *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(MGparm));
    VASSERT( thee != VNULL);
    VASSERT( MGparm_ctor2(thee, type) );

    return thee;
}

VPUBLIC int MGparm_ctor2(MGparm *thee, MGparm_CalcType type) {

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
    thee->setchgm = 0;

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
    thee->setrank = 0;
    thee->setsize = 0;
    thee->setofrac = 0;
    for (i=0; i<6; i++) thee->partDisjOwnSide[i] = 0;
	thee->setasync = 0;

    /* *** Default parameters for TINKER *** */
    thee->chgs = VCM_CHARGE;

    return 1; 
}

VPUBLIC void MGparm_dtor(MGparm **thee) {
    if ((*thee) != VNULL) {
        MGparm_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(MGparm), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void MGparm_dtor2(MGparm *thee) { ; }

VPUBLIC int MGparm_check(MGparm *thee) { 

    int rc, i, tdime[3], ti, tnlev[3], nlev;

    rc = VRC_SUCCESS;
	
	Vnm_print(0, "MGparm_check:  checking MGparm object of type %d.\n", 
			  thee->type);

    /* Check to see if we were even filled... */
    if (!thee->parsed) {
        Vnm_print(2, "MGparm_check:  not filled!\n");
        return VRC_FAILURE;
    }

    /* Check generic settings */
    if (!thee->setdime) {
        Vnm_print(2, "MGparm_check:  DIME not set!\n");
        rc = VRC_FAILURE;
    }
    if (!thee->setchgm) {
        Vnm_print(2, "MGparm_check: CHGM not set!\n");
        return VRC_FAILURE;
    }


    /* Check sequential manual & dummy settings */
    if ((thee->type == MCT_MANUAL) || (thee->type == MCT_DUMMY)) {
        if ((!thee->setgrid) && (!thee->setglen)) {
            Vnm_print(2, "MGparm_check:  Neither GRID nor GLEN set!\n");
            rc = VRC_FAILURE;
        }
        if ((thee->setgrid) && (thee->setglen)) {
            Vnm_print(2, "MGparm_check:  Both GRID and GLEN set!\n");
            rc = VRC_FAILURE;
        }
        if (!thee->setgcent) {
            Vnm_print(2, "MGparm_check:  GCENT not set!\n");
            rc = VRC_FAILURE;
        }
    }
 
    /* Check sequential and parallel automatic focusing settings */
    if ((thee->type == MCT_AUTO) || (thee->type == MCT_PARALLEL)) {
        if (!thee->setcglen) {
            Vnm_print(2, "MGparm_check:  CGLEN not set!\n");
            rc = VRC_FAILURE;
        }
        if (!thee->setfglen) {
            Vnm_print(2, "MGparm_check:  FGLEN not set!\n");
            rc = VRC_FAILURE;
        }
        if (!thee->setcgcent) {
            Vnm_print(2, "MGparm_check:  CGCENT not set!\n");
            rc = VRC_FAILURE;
        }
        if (!thee->setfgcent) {
            Vnm_print(2, "MGparm_check:  FGCENT not set!\n");
            rc = VRC_FAILURE;
        }
    }

    /* Check parallel automatic focusing settings */
    if (thee->type == MCT_PARALLEL) {
        if (!thee->setpdime) {
            Vnm_print(2, "MGparm_check:  PDIME not set!\n");
            rc = VRC_FAILURE;
        }
        if (!thee->setrank) {
            Vnm_print(2, "MGparm_check:  PROC_RANK not set!\n");
            rc = VRC_FAILURE;
        }
        if (!thee->setsize) {
            Vnm_print(2, "MGparm_check:  PROC_SIZE not set!\n");
            rc = VRC_FAILURE;
        }
        if (!thee->setofrac) {
            Vnm_print(2, "MGparm_check:  OFRAC not set!\n");
            rc = VRC_FAILURE;
        }
    }
 
    /* Perform a sanity check on nlev and dime, resetting values as necessary */
    if (rc == 1) {
	/* Calculate the actual number of grid points and nlev to satisfy the
	 * formula:  n = c * 2^(l+1) + 1, where n is the number of grid points,
	 * c is an integer, and l is the number of levels */
        for (i=0; i<3; i++) {
            /* See if the user picked a reasonable value, if not then fix it */
            ti = thee->dime[i] - 1;
            if (ti == VPOW(2, (thee->nlev+1))) {
                tnlev[i] = thee->nlev;
                tdime[i] = thee->dime[i];
            } else {
                tdime[i] = thee->dime[i];
                ti = tdime[i] - 1;
                tnlev[i] = 0;
                /* Find the maximum number of times this dimension can be
                 * divided by two */
                while (VEVEN(ti)) {
                    (tnlev[i])++;
                    ti = (int)ceil(0.5*ti);
                }
                (tnlev[i])--;
                /* We'd like to have at least VMGNLEV levels in the multigrid
                 * hierarchy.  This means that the dimension needs to be
                 * c*2^VMGNLEV + 1, where c is an integer. */
                if ((tdime[i] > 65) && (tnlev[i] < VMGNLEV)) {
                    Vnm_print(2, "NOsh:  Bad dime[%d]  = %d (%d nlev)!\n",
                      i, tdime[i], tnlev[i]);
                    ti = (int)(tdime[i]/VPOW(2.,(VMGNLEV+1)));
                    if (ti < 1) ti = 1;
                    tdime[i] = ti*(int)(VPOW(2.,(VMGNLEV+1))) + 1;
                    tnlev[i] = 4;
                    Vnm_print(2, "NOsh:  Reset dime[%d] to %d and (nlev = %d).\n", i, tdime[i], VMGNLEV);
                }
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
		if (thee->nlev <= 0) {
			Vnm_print(2, "MGparm_check:  illegal nlev (%d); check your grid dimensions!\n", thee->nlev);
			rc = VRC_FAILURE;
		}
		if (thee->nlev < 2) {
			Vnm_print(2, "MGparm_check:  you're using a very small nlev (%d) and therefore\n", thee->nlev);
			Vnm_print(2, "MGparm_check:  will not get the optimal performance of the multigrid\n");
			Vnm_print(2, "MGparm_check:  algorithm.  Please check your grid dimensions.\n");
		}
        for (i=0; i<3; i++) thee->dime[i] = tdime[i];
    }

    return rc;
}

VPUBLIC void MGparm_copy(MGparm *thee, MGparm *parm) {

    int i;

    VASSERT(thee != VNULL);
    VASSERT(parm != VNULL);

    
    thee->type = parm->type;
    thee->parsed = parm->parsed;

    /* *** GENERIC PARAMETERS *** */
    for (i=0; i<3; i++) thee->dime[i] = parm->dime[i];
    thee->setdime = parm->setdime;
    thee->chgm = parm->chgm;
    thee->setchgm = parm->setchgm;
    thee->chgs = parm->chgs;

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
      thee->partDisjCenter[i] = parm->partDisjCenter[i];
    for (i=0; i<3; i++) 
      thee->partDisjLength[i] = parm->partDisjLength[i];
    for (i=0; i<6; i++) 
      thee->partDisjOwnSide[i] = parm->partDisjOwnSide[i];
    for (i=0; i<3; i++) thee->pdime[i] = parm->pdime[i];
    thee->setpdime = parm->setpdime;
    thee->proc_rank = parm->proc_rank;
    thee->setrank = parm->setrank;
    thee->proc_size = parm->proc_size;
    thee->setsize = parm->setsize;
    thee->ofrac = parm->ofrac;
    thee->setofrac = parm->setofrac;
	thee->setasync = parm->setasync;
	thee->async = parm->async;

}

VPRIVATE int MGparm_parseDIME(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    int ti;

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

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int MGparm_parseCHGM(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    int ti;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%d", &ti) == 1) {
        thee->chgm = ti;
        thee->setchgm = 1;
        Vnm_print(2, "NOsh:  Warning -- parsed deprecated statment \"chgm %d\".\n", ti);
        Vnm_print(2, "NOsh:  Please use \"chgm ");
        switch (thee->chgm) {
            case VCM_TRIL:
                Vnm_print(2, "spl0");
                break;
            case VCM_BSPL2:
                Vnm_print(2, "spl2");
                break;
            case VCM_BSPL4:
                Vnm_print(2, "spl4"); 
                break;
            default:
                Vnm_print(2, "UNKNOWN");
                break;
        }
        Vnm_print(2, "\" instead!\n");
        return 1;
    } else if (Vstring_strcasecmp(tok, "spl0") == 0) {
        thee->chgm = VCM_TRIL;
        thee->setchgm = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "spl2") == 0) {
        thee->chgm = VCM_BSPL2;
        thee->setchgm = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "spl4") == 0) {
        thee->chgm = VCM_BSPL4;
        thee->setchgm = 1;
        return 1;
    } else {
        Vnm_print(2, "NOsh:  Unrecognized parameter (%s) when parsing \
chgm!\n", tok);
        return -1;
    }
    return -1;

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int MGparm_parseNLEV(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    int ti;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%d", &ti) == 0) {
        Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing NLEV \
keyword!\n", tok);
        return -1;
    } else thee->nlev = ti;
    thee->setnlev = 1;
    return 1;

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int MGparm_parseGRID(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    double tf;

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

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int MGparm_parseGLEN(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    double tf;

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

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int MGparm_parseGAMMA(MGparm *thee, Vio *sock) {
	
    char tok[VMAX_BUFSIZE];
	
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	Vnm_print(2, "parseMG:  GAMMA keyword deprecated!\n");
	Vnm_print(2, "parseMG:  Please see new APOLAR documentation.\n");
	return 1;
	
VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
	return -1;
}

VPRIVATE int MGparm_parseGCENT(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    double tf;
    int ti;

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
                thee->cmeth = MCM_MOLECULE;
				/* Subtract 1 here to convert user numbering (1, 2, 3, ...) into
				array index */
                thee->centmol = ti - 1;
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

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int MGparm_parseCGLEN(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    double tf;

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

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int MGparm_parseFGLEN(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    double tf;

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

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int MGparm_parseCGCENT(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    double tf;
    int ti;

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
				thee->ccmeth = MCM_MOLECULE;
				/* Subtract 1 here to convert user numbering (1, 2, 3, ...) into 
				array index */
                thee->ccentmol = ti - 1;
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

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int MGparm_parseFGCENT(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    double tf;
    int ti;

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
                thee->fcmeth = MCM_MOLECULE;
				/* Subtract 1 here to convert user numbering (1, 2, 3, ...) into
				array index */
                thee->fcentmol = ti - 1;
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

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int MGparm_parsePDIME(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    int ti;

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

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int MGparm_parseOFRAC(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-int (%s) while parsing OFRAC \
keyword!\n", tok);
        return -1;
    }
    thee->ofrac = tf;
    thee->setofrac = 1;
    return 1;

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int MGparm_parseASYNC(MGparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    int ti;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%i", &ti) == 0) {
        Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing ASYNC \
keyword!\n", tok);
        return -1;
    }
    thee->async = ti;
    thee->setasync = 1;
    return 1;

    VERROR1:
        Vnm_print(2, "parseMG:  ran out of tokens!\n");
        return -1;
}

VPUBLIC int MGparm_parseToken(MGparm *thee, char tok[VMAX_BUFSIZE], 
  Vio *sock) {

    if (thee == VNULL) {
        Vnm_print(2, "parseMG:  got NULL thee!\n"); 
        return -1;
    }
    if (sock == VNULL) {
        Vnm_print(2, "parseMG:  got NULL socket!\n"); 
        return -1;
    }

    if (Vstring_strcasecmp(tok, "dime") == 0) {
        return MGparm_parseDIME(thee, sock);
    } else if (Vstring_strcasecmp(tok, "chgm") == 0) {
        return MGparm_parseCHGM(thee, sock);
    } else if (Vstring_strcasecmp(tok, "nlev") == 0) {
        return MGparm_parseNLEV(thee, sock);
    } else if (Vstring_strcasecmp(tok, "grid") == 0) {
        return MGparm_parseGRID(thee, sock);
    } else if (Vstring_strcasecmp(tok, "glen") == 0) {
        return MGparm_parseGLEN(thee, sock);
    } else if (Vstring_strcasecmp(tok, "gcent") == 0) {
        return MGparm_parseGCENT(thee, sock);
    } else if (Vstring_strcasecmp(tok, "cglen") == 0) {
        return MGparm_parseCGLEN(thee, sock);
    } else if (Vstring_strcasecmp(tok, "fglen") == 0) {
        return MGparm_parseFGLEN(thee, sock);
    } else if (Vstring_strcasecmp(tok, "cgcent") == 0) {
        return MGparm_parseCGCENT(thee, sock);
    } else if (Vstring_strcasecmp(tok, "fgcent") == 0) {
        return MGparm_parseFGCENT(thee, sock);
    } else if (Vstring_strcasecmp(tok, "pdime") == 0) {
        return MGparm_parsePDIME(thee, sock);
    } else if (Vstring_strcasecmp(tok, "ofrac") == 0) {
        return MGparm_parseOFRAC(thee, sock);
    } else if (Vstring_strcasecmp(tok, "async") == 0) {
        return MGparm_parseASYNC(thee, sock);
	} else if (Vstring_strcasecmp(tok, "gamma") == 0) {
        return MGparm_parseGAMMA(thee, sock);
    } else {
        Vnm_print(2, "parseMG:  Unrecognized keyword (%s)!\n", tok);
        return -1;
    }
   
    return 0;

}
