/**
 *  @file    pbeparm.c
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @brief   Class PBEparm methods
 *  @version $Id$
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
 * Copyright (c) 1999-2002. The Regents of the University of California
 *                          (Regents).  All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above copyright notice, this paragraph and the
 * following two paragraphs appear in all copies, modifications, and
 * distributions.
 *
 * IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
 * REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
 * ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
 * TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS.
 *
 * @endverbatim
 */

#include "apbscfg.h"
#include "apbs/pbeparm.h"
#include "apbs/vstring.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class PBEparm: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_MGPARM)

#endif /* if !defined(VINLINE_MGPARM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class PBEparm: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  PBEparm_ctor
//
// Author: Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC PBEparm* PBEparm_ctor() {

    /* Set up the structure */
    PBEparm *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(PBEparm));
    VASSERT( thee != VNULL);
    VASSERT( PBEparm_ctor2(thee) );

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  PBEparm_ctor2
//
// Purpose:  Construct the PBEparm object
//
// Notes:    Constructor broken into two parts for FORTRAN users.
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int PBEparm_ctor2(PBEparm *thee) {

    int i;

    if (thee == VNULL) return 0;

    thee->parsed = 0;

    thee->setmolid = 0;
    thee->setnonlin = 0;
    thee->setbcfl = 0;
    thee->setnion = 0;
    for (i=0; i<MAXION; i++) thee->setion[i] = 0;
    thee->setpdie = 0;
    thee->setsdie = 0;
    thee->setsrfm = 0;
    thee->setsrad = 0;
    thee->setswin = 0; 
    thee->settemp = 0;
    thee->setgamma = 0;
    thee->setcalcenergy = 0;      
    thee->setcalcforce = 0;       
    thee->setwritepot = 0; 
    thee->setwriteacc = 0; 
    thee->nion = 0;
    thee->swin = 0;
    thee->srad = 1.4;

    return 1; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  PBEparm_dtor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void PBEparm_dtor(PBEparm **thee) {
    if ((*thee) != VNULL) {
        PBEparm_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(PBEparm), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  PBEparm_dtor2
//
// Purpose:  Destroy the atom object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void PBEparm_dtor2(PBEparm *thee) { ; }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  PBEparm_check
//
// Purpose:  Check the parameter settings for internal consistency
//
// Returns:  1 if OK, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int PBEparm_check(PBEparm *thee) { 

    int i;

    /* Check to see if we were even filled... */
    if (!thee->parsed) {
        Vnm_print(2, "PBEparm_check:  not filled!\n");
        return 0;
    }

    if (!thee->setmolid) {
        Vnm_print(2, "PBEparm_check:  MOL not set!\n");
        return 0;
    }
    if (!thee->setnonlin) {
        Vnm_print(2, "PBEparm_check:  LPBE or NPBE not set!\n");
        return 0;
    }
    if (!thee->setbcfl) {
        Vnm_print(2, "PBEparm_check:  BCFL not set!\n");
        return 0;
    }
    if (!thee->setnion) {
        thee->setnion = 1;
        thee->nion = 0;
    } 
    for (i=0; i<thee->nion; i++) {
        if (!thee->setion[i]) {
            Vnm_print(2, "PBEparm_check: ION #%d not set!\n",i);
            return 0;
        }
    }
    if (!thee->setpdie) {
        Vnm_print(2, "PBEparm_check: PDIE not set!\n");
        return 0;
    }
    if (!thee->setsdie) {
        Vnm_print(2, "PBEparm_check: SDIE not set!\n");
        return 0;
    }
    if (!thee->setsrfm) {
        Vnm_print(2, "PBEparm_check: SRFM not set!\n");
        return 0;
    }
    if (((thee->srfm==0) || (thee->srfm==1)) && (!thee->setsrad)) {
        Vnm_print(2, "PBEparm_check: SRAD not set!\n");
        return 0;
    }
    if ((thee->srfm==2) && (!thee->setswin)) {
        Vnm_print(2, "PBEparm_check: SWIN not set!\n");
        return 0;
    }
    if (!thee->settemp) {
        Vnm_print(2, "PBEparm_check: TEMP not set!\n");
        return 0;
    }
    if (!thee->setgamma) {
        Vnm_print(2, "PBEparm_check: GAMMA not set!\n");
        return 0;
    }
    if (!thee->setcalcenergy) thee->calcenergy = 0;
    if (!thee->setcalcforce) thee->calcforce = 0;
    if (!thee->setwritepot) thee->writepot = 0;
    if (!thee->setwriteacc) thee->writeacc = 0;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  PBEparm_copy
//  
// Purpose:  Copy parm into thee
//
// Args:     parm    object to copy into thee
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void PBEparm_copy(PBEparm *thee, PBEparm *parm) {

    int i;

    VASSERT(thee != VNULL);
    VASSERT(parm != VNULL);

    thee->molid = parm->molid;
    thee->setmolid = parm->setmolid;
    thee->nonlin = parm->nonlin; 
    thee->setnonlin = parm->setnonlin;
    thee->bcfl = parm->bcfl;
    thee->setbcfl = parm->setbcfl;
    thee->nion = parm->nion;
    thee->setnion = parm->setnion;
    for (i=0; i<MAXION; i++) {
        thee->ionq[i] = parm->ionq[i];
        thee->ionc[i] = parm->ionc[i];
        thee->ionr[i] = parm->ionr[i];
        thee->setion[i] = parm->setion[i];
    };
    thee->pdie = parm->pdie;
    thee->setpdie = parm->setpdie;
    thee->sdie = parm->sdie;
    thee->setsdie = parm->setsdie;
    thee->srfm = parm->srfm;
    thee->setsrfm = parm->setsrfm;
    thee->srad = parm->srad;
    thee->setsrad = parm->setsrad;
    thee->swin = parm->swin;
    thee->setswin = parm->setswin;
    thee->temp = parm->temp;
    thee->settemp = parm->settemp;
    thee->gamma = parm->gamma;
    thee->setgamma = parm->setgamma;
    thee->calcenergy = parm->calcenergy;
    thee->setcalcenergy = parm->setcalcenergy;
    thee->calcforce = parm->calcforce;
    thee->setcalcforce = parm->setcalcforce;
    thee->writepot = parm->writepot;
    thee->setwritepot = parm->setwritepot;
    for (i=0; i<VMAX_ARGLEN; i++) thee->writepotstem[i] = parm->writepotstem[i];
    thee->writepotfmt = parm->writepotfmt;
    thee->writeacc = parm->writeacc;
    thee->setwriteacc = parm->setwriteacc;
    for (i=0; i<VMAX_ARGLEN; i++) thee->writeaccstem[i] = parm->writeaccstem[i];
    thee->writeaccfmt = parm->writeaccfmt;
    thee->writemat = parm->writemat;
    thee->setwritemat = parm->setwritemat;
    for (i=0; i<VMAX_ARGLEN; i++) thee->writematstem[i] = parm->writematstem[i];
    thee->writematflag = parm->writematflag;

    thee->parsed = parm->parsed;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  PBEparm_parseToken
//
// Purpose:  Parse a PBE keyword
//
// Returns:  1 if matched and assigned
//          -1 if matched, but there's some sort of error (i.e., too few args)
//           0 if not matched
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int PBEparm_parseToken(PBEparm *thee, char tok[VMAX_BUFSIZE], 
  Vio *sock) {

    int ti;
    double tf;

    if (thee == VNULL) {
        Vnm_print(2, "parsePBE:  got NULL thee!\n");
        return -1;
    }
    if (sock == VNULL) {
        Vnm_print(2, "parsePBE:  got NULL socket!\n");
        return -1;
    }

    if (Vstring_strcasecmp(tok, "mol") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing MOL \
keyword!\n", tok);
            return -1;
        } 
        thee->molid = ti;
        thee->setmolid = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "lpbe") == 0) {
        thee->nonlin = 0;
        thee->setnonlin = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "npbe") == 0) {
        thee->nonlin = 1;
        thee->setnonlin = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "bcfl") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing BCFL \
keyword!\n", tok);
            return -1;
        } 
        thee->bcfl = ti;
        thee->setbcfl = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "ion") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION \
keyword!\n", tok);
            return -1;
        }
        thee->ionq[thee->nion] = tf;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION \
keyword!\n", tok);
            return -1;
        }
        thee->ionc[thee->nion] = tf;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION \
keyword!\n", tok);
            return -1;
        }
        thee->ionr[thee->nion] = tf;
        thee->setion[thee->nion] = 1;
        (thee->nion)++;
        thee->setnion = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "pdie") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing PDIE \
keyword!\n", tok);
            return -1;
        }
        thee->pdie = tf;
        thee->setpdie = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "sdie") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SDIE \
keyword!\n", tok);
            return -1;
        }
        thee->sdie = tf;
        thee->setsdie = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "srfm") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing SRFM \
keyword!\n", tok);
            return -1;
        }
        thee->srfm = ti;
        thee->setsrfm = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "srad") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SRAD \
keyword!\n", tok);
            return -1;
        }
        thee->srad = tf;
        thee->setsrad = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "swin") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
           Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SWIN \
keyword!\n", tok);
           return -1;
        }
        thee->swin = tf;
        thee->setswin = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "temp") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing TEMP \
keyword!\n", tok);
            return -1;
        }
        thee->temp = tf;
        thee->settemp = 1; 
        return 1;
    } else if (Vstring_strcasecmp(tok, "gamma") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GAMMA \
keyword!\n", tok);
            return -1;
        }
        thee->gamma = tf;
        thee->setgamma = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "calcenergy") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing \
WRITEENERGY keyword!\n", tok);
            return -1;
        }
        thee->calcenergy = ti;
        thee->setcalcenergy = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "calcforce") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing \
WRITEFORCE keyword!\n", tok);
            return -1;
        }
        thee->calcforce = ti;
        thee->setcalcforce = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "writepot") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing WRITEPOT \
keyword!\n", tok);
            return -1;
        }
        thee->writepot = ti;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (Vstring_strcasecmp(tok, "dx") == 0) {
            thee->writepotfmt = 0;
        } else if (Vstring_strcasecmp(tok, "avs") == 0) {
            thee->writepotfmt = 1;
        } else if (Vstring_strcasecmp(tok, "uhbd") == 0) {
            thee->writepotfmt = 2;
        } else {
            Vnm_print(2, "NOsh:  Invalid format (%s) while parsing \
WRITEPOT keyword!\n", tok);
            return -1;
        }
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        strncpy(thee->writepotstem, tok, VMAX_ARGLEN);
        thee->setwritepot = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "writeacc") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing WRITEACC \
keyword!\n", tok);
            return -1;
        } 
        thee->writeacc = ti;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (Vstring_strcasecmp(tok, "dx") == 0) {
            thee->writeaccfmt = 0;
        } else if (Vstring_strcasecmp(tok, "avs") == 0) {
            thee->writeaccfmt = 1;
        } else if (Vstring_strcasecmp(tok, "uhbd") == 0) {
            thee->writeaccfmt = 2;
        } else {
            Vnm_print(2, "NOsh:  Invalid format (%s) while parsing \
WRITEACC keyword!\n", tok);
            return -1;
        }
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        strncpy(thee->writeaccstem, tok, VMAX_ARGLEN);
        thee->setwriteacc = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "writemat") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing WRITEMAT \
keyword!\n", tok);
            return -1;
        }
        thee->writemat = ti;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (Vstring_strcasecmp(tok, "poisson") == 0) {
            thee->writematflag = 0;
        } else if (Vstring_strcasecmp(tok, "full") == 0) {
            thee->writematflag = 1;
        } else {
            Vnm_print(2, "NOsh:  Invalid format (%s) while parsing \
WRITEMAT keyword!\n", tok);
            return -1;
        }
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        strncpy(thee->writematstem, tok, VMAX_ARGLEN);
        thee->setwriteacc = 1;
        return 1;

    }

    return 0;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;

}
