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
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Center for Computational Biology
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2003.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2003.  The Regents of the University of
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

VPUBLIC double PBEparm_getIonCharge(PBEparm *thee, int i) {
    VASSERT(thee != VNULL);
	VASSERT(i < thee->nion);
    return thee->ionq[i];
}
VPUBLIC double PBEparm_getIonConc(PBEparm *thee, int i) {
    VASSERT(thee != VNULL);
    VASSERT(i < thee->nion);
    return thee->ionc[i];
}
VPUBLIC double PBEparm_getIonRadius(PBEparm *thee, int i) {
    VASSERT(thee != VNULL);
    VASSERT(i < thee->nion);
    return thee->ionr[i];
}


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
    thee->numwrite = 0; 
    thee->setwritemat = 0; 
    thee->nion = 0;
    thee->swin = 0;
    thee->srad = 1.4;
    thee->useDielMap = 0;
    thee->useKappaMap = 0;
    thee->useChargeMap = 0;

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
    if (((thee->srfm==VSM_MOL) || (thee->srfm==VSM_MOLSMOOTH)) \
      && (!thee->setsrad)) {
        Vnm_print(2, "PBEparm_check: SRAD not set!\n");
        return 0;
    }
    if ((thee->srfm==VSM_SPLINE) && (!thee->setswin)) {
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
    if (!thee->setwritemat) thee->writemat = 0;

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

    int i, j;

    VASSERT(thee != VNULL);
    VASSERT(parm != VNULL);

    thee->molid = parm->molid;
    thee->setmolid = parm->setmolid;
    thee->useDielMap = parm->useDielMap;
    thee->dielMapID = parm->dielMapID;
    thee->useKappaMap = parm->useKappaMap;
    thee->kappaMapID = parm->kappaMapID;
    thee->useChargeMap = parm->useChargeMap;
    thee->chargeMapID = parm->chargeMapID;
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
    thee->numwrite = parm->numwrite;
    for (i=0; i<PBEPARM_MAXWRITE; i++) {
        thee->writetype[i] = parm->writetype[i];
        thee->writefmt[i] = parm->writefmt[i];
        for (j=0; j<VMAX_ARGLEN; j++) 
          thee->writestem[i][j] = parm->writestem[i][j];
    }
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
    Vdata_Type writetype;
    Vdata_Format writefmt;

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
        /* We can either parse int flag... */
        if (sscanf(tok, "%d", &ti) == 1) {
            thee->bcfl = ti;
            thee->setbcfl = 1;
            return 1;
        /* ...or the word */
        } else {
            if (Vstring_strcasecmp(tok, "zero") == 0) {
                thee->bcfl = BCFL_ZERO;
                thee->setbcfl = 1;
                return 1;
            } else if (Vstring_strcasecmp(tok, "sdh") == 0) {
                thee->bcfl = BCFL_SDH;
                thee->setbcfl = 1;
                return 1;
            } else if (Vstring_strcasecmp(tok, "mdh") == 0) {
                thee->bcfl = BCFL_MDH;
                thee->setbcfl = 1;
                return 1;
            } else if (Vstring_strcasecmp(tok, "focus") == 0) {
                thee->bcfl = BCFL_FOCUS;
                thee->setbcfl = 1;
                return 1;
            } else {
                Vnm_print(2, "NOsh:  parsed unknown BCFL parameter (%s)!\n",
                  tok);
                return -1;
            }
        }
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
        if (sscanf(tok, "%d", &ti) == 1) {
            thee->srfm = ti;
            thee->setsrfm = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "mol") == 0) {
            thee->srfm = VSM_MOL;
            thee->setsrfm = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "smol") == 0) {
            thee->srfm = VSM_MOLSMOOTH;
            thee->setsrfm = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "spl2") == 0) {
            thee->srfm = VSM_MOLSMOOTH;
            thee->setsrfm = 1;
            return 1;
        } else {
            Vnm_print(2, "NOsh:  Unrecongnized keyword (%s) when parsing \
srfm!\n", tok);
            return -1;
        }
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
    } else if (Vstring_strcasecmp(tok, "usemap") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        Vnm_print(0, "PBEparm_parseToken:  Read %s...\n", tok);
        if (Vstring_strcasecmp(tok, "diel") == 0) {
            thee->useDielMap = 1;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1); 
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-int (%s) while parsing \
USEMAP DIEL keyword!\n", tok);
                return -1;
            } 
            thee->dielMapID = ti;
            return 1;
        } else if (Vstring_strcasecmp(tok, "kappa") == 0) {
            thee->useKappaMap = 1;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-int (%s) while parsing \
USEMAP KAPPA keyword!\n", tok);
                return -1;
            }
            thee->kappaMapID = ti;
            return 1;
        } else if (Vstring_strcasecmp(tok, "charge") == 0) {
            thee->useChargeMap = 1;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) { 
                Vnm_print(2, "NOsh:  Read non-int (%s) while parsing \
USEMAP CHARGE keyword!\n", tok);
                return -1;
            }
            thee->chargeMapID = ti;
            return 1;
        } else {
            Vnm_print(2, "NOsh:  Read undefined keyword (%s) while parsing \
USEMAP statement!\n", tok);
            return -1;
        }
    } else if (Vstring_strcasecmp(tok, "calcenergy") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 1) {
            return -1;
            thee->calcenergy = ti;
            thee->setcalcenergy = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "no") == 0) {
            thee->calcenergy = 0;
            thee->setcalcenergy = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "total") == 0) {
            thee->calcenergy = 1;
            thee->setcalcenergy = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "comps") == 0) {
            thee->calcenergy = 2;
            thee->setcalcenergy = 1;
            return 1;
        } else {
            Vnm_print(2, "NOsh:  Unrecognized parameter (%s) while parsing \
calcenergy!\n", tok);
            return -1;
        }
    } else if (Vstring_strcasecmp(tok, "calcforce") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 1) {
            thee->calcforce = ti;
            thee->setcalcforce = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "no") == 0) {
            thee->calcforce = 0;
            thee->setcalcforce = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "total") == 0) {
            thee->calcforce = 1;
            thee->setcalcforce = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "comps") == 0) {
            thee->calcforce = 2;
            thee->setcalcforce = 1;
            return 1;
        } else {
            Vnm_print(2, "NOsh:  Unrecognized parameter (%s) while parsing \
calcforce!\n", tok);
            return -1;
        }
    } else if (Vstring_strcasecmp(tok, "write") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (Vstring_strcasecmp(tok, "pot") == 0) {
            writetype = VDT_POT;
        } else if (Vstring_strcasecmp(tok, "charge") == 0) {
            writetype = VDT_CHARGE;
        } else if (Vstring_strcasecmp(tok, "smol") == 0) {
            writetype = VDT_SMOL;
        } else if (Vstring_strcasecmp(tok, "dielx") == 0) {
            writetype = VDT_DIELX;
        } else if (Vstring_strcasecmp(tok, "diely") == 0) {
            writetype = VDT_DIELY;
        } else if (Vstring_strcasecmp(tok, "dielz") == 0) {
            writetype = VDT_DIELZ;
        } else if (Vstring_strcasecmp(tok, "kappa") == 0) {
            writetype = VDT_KAPPA;
        } else if (Vstring_strcasecmp(tok, "sspl") == 0) {
            writetype = VDT_SSPL;
        } else if (Vstring_strcasecmp(tok, "vdw") == 0) {
            writetype = VDT_VDW;
        } else if (Vstring_strcasecmp(tok, "ivdw") == 0) {
            writetype = VDT_IVDW;
        } else if (Vstring_strcasecmp(tok, "lap") == 0) {
            writetype = VDT_LAP;
        } else if (Vstring_strcasecmp(tok, "edens") == 0) {
            writetype = VDT_EDENS;
        } else if (Vstring_strcasecmp(tok, "ndens") == 0) {
            writetype = VDT_NDENS;
        } else if (Vstring_strcasecmp(tok, "qdens") == 0) {
            writetype = VDT_QDENS;
        } else {
            Vnm_print(2, "PBEparm_parse:  Invalid data type (%s) to write!\n",
               tok);
            return -1;
        }
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (Vstring_strcasecmp(tok, "dx") == 0) {
            writefmt = VDF_DX;
        } else if (Vstring_strcasecmp(tok, "uhbd") == 0) {
            writefmt = VDF_UHBD;
        } else if (Vstring_strcasecmp(tok, "avs") == 0) {
            writefmt = VDF_AVS;
        } else {
            Vnm_print(2, "PBEparm_parse:  Invalid data format (%s) to write!\n",
               tok);
            return -1;
        }
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        strncpy(thee->writestem[thee->numwrite], tok, VMAX_ARGLEN);
        thee->writetype[thee->numwrite] = writetype;
        thee->writefmt[thee->numwrite] = writefmt;
        (thee->numwrite)++;
        return 1;

    } else if (Vstring_strcasecmp(tok, "writemat") == 0) {
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
        thee->setwritemat = 1;
        thee->writemat = 1;
        return 1;

    }

    return 0;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;

}
