/**
 *  @file    femparm.c
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @brief   Class FEMparm methods
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
#include "apbs/femparm.h"

#if !defined(VINLINE_MGPARM)

#endif /* if !defined(VINLINE_MGPARM) */

VPUBLIC FEMparm* FEMparm_ctor(FEMparm_CalcType type) {

    /* Set up the structure */
    FEMparm *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(FEMparm));
    VASSERT( thee != VNULL);
    VASSERT( FEMparm_ctor2(thee, type) );

    return thee;
}

VPUBLIC int FEMparm_ctor2(FEMparm *thee, FEMparm_CalcType type) {

    if (thee == VNULL) return 0;

    thee->parsed = 0;
    thee->type = type;
    thee->settype = 1;

    thee->setdomainLength = 0;
    thee->setetol = 0;
    thee->setekey = 0;
    thee->setakeyPRE = 0;
    thee->setakeySOLVE = 0;
    thee->settargetNum = 0;
    thee->settargetRes = 0;
    thee->setmaxsolve = 0;
    thee->setmaxvert = 0;

    return 1; 
}

VPUBLIC void FEMparm_dtor(FEMparm **thee) {
    if ((*thee) != VNULL) {
        FEMparm_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(FEMparm), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void FEMparm_dtor2(FEMparm *thee) { ; }

VPUBLIC int FEMparm_check(FEMparm *thee) { 

    int rc;
    rc = 1;

    if (!thee->parsed) {
        Vnm_print(2, "FEMparm_check:  not filled!\n");
        return 0;
    }
    if (!thee->settype) {
        Vnm_print(2, "FEMparm_check:  type not set!\n");
        rc = 0;
    }
    if (!thee->setdomainLength) {
        Vnm_print(2, "FEMparm_check:  domainLength not set!\n");
        rc = 0;
    }
    if (!thee->setetol) {
        Vnm_print(2, "FEMparm_check:  etol not set!\n");
        rc = 0;
    }
    if (!thee->setekey) {
        Vnm_print(2, "FEMparm_check:  ekey not set!\n");
        rc = 0;
    }
    if (!thee->setakeyPRE) {
        Vnm_print(2, "FEMparm_check:  akeyPRE not set!\n");
        rc = 0;
    }
    if (!thee->setakeySOLVE) {
        Vnm_print(2, "FEMparm_check:  akeySOLVE not set!\n");
        rc = 0;
    }
    if (!thee->settargetNum) {
        Vnm_print(2, "FEMparm_check:  targetNum not set!\n");
        rc = 0;
    }
    if (!thee->settargetRes) {
        Vnm_print(2, "FEMparm_check:  targetRes not set!\n");
        rc = 0;
    }
    if (!thee->setmaxsolve) {
        Vnm_print(2, "FEMparm_check:  maxsolve not set!\n");
        rc = 0;
    }
    if (!thee->setmaxvert) {
        Vnm_print(2, "FEMparm_check:  maxvert not set!\n");
        rc = 0;
    }

    return rc;
}

VPUBLIC int FEMparm_parseToken(FEMparm *thee, char tok[VMAX_BUFSIZE],
  Vio *sock) {

    int i, ti;
    double tf;

    if (thee == VNULL) {
        Vnm_print(2, "parseFE:  got NULL thee!\n");
        return -1;
    }

    if (sock == VNULL) {
        Vnm_print(2, "parseFE:  got NULL socket!\n");
        return -1;
    }

    if (Vstring_strcasecmp(tok, "domainLength") == 0) {
        for (i=0; i<3; i++) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "parseFE:  Read non-double (%s) while parsing \
DOMAINRADIUS keyword!\n", tok);
                return -1;
            } 
            thee->domainLength[i] = tf;
        }
        thee->setdomainLength = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "etol") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "parseFE:  Read non-double (%s) while parsing \
ETOL keyword!\n", tok);
            return -1;
        }
        thee->etol = tf;
        thee->setetol = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "ekey") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (Vstring_strcasecmp(tok, "simp") == 0) {
            thee->ekey = FET_SIMP;
            thee->setekey = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "glob") == 0) {
            thee->ekey = FET_GLOB;
            thee->setekey = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "frac") == 0) {
            thee->ekey = FET_FRAC;
            thee->setekey = 1;
            return 1;
        } else {
            Vnm_print(2, "parseFE:  undefined value (%s) for ekey!\n", tok);
            return -1;
        }
    } else if (Vstring_strcasecmp(tok, "akeyPRE") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (Vstring_strcasecmp(tok, "unif") == 0) {
            thee->akeyPRE = FRT_UNIF;
            thee->setakeyPRE = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "geom") == 0) {
            thee->akeyPRE = FRT_GEOM;
            thee->setakeyPRE = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "resi") == 0) {
            thee->akeyPRE = FRT_RESI;
            thee->setakeyPRE = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "dual") == 0) {
            thee->akeyPRE = FRT_DUAL;
            thee->setakeyPRE = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "loca") == 0) {
            thee->akeyPRE = FRT_LOCA;
            thee->setakeyPRE = 1;
            return 1;
        } else {
            Vnm_print(2, "parseFE:  undefined value (%s) for akeyPRE!\n", tok);
            return -1;
        }
    } else if (Vstring_strcasecmp(tok, "akeySOLVE") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (Vstring_strcasecmp(tok, "unif") == 0) {
            thee->akeySOLVE = FRT_UNIF;
            thee->setakeySOLVE = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "geom") == 0) {
            thee->akeySOLVE = FRT_GEOM;
            thee->setakeySOLVE = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "resi") == 0) {
            thee->akeySOLVE = FRT_RESI;
            thee->setakeySOLVE = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "dual") == 0) {
            thee->akeySOLVE = FRT_DUAL;
            thee->setakeySOLVE = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "loca") == 0) {
            thee->akeySOLVE = FRT_LOCA;
            thee->setakeySOLVE = 1;
            return 1;
        } else {
            Vnm_print(2, "parseFE:  undefined value (%s) for akeyPRE!\n", tok);
            return -1;
        }
    } else if (Vstring_strcasecmp(tok, "targetNum") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "parseFE:  read non-int (%s) for targetNum!\n", tok);
            return -1;
        }
        thee->targetNum = ti;
        thee->settargetNum = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "targetRes") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%f", &tf) == 0) {
            Vnm_print(2, "parseFE:  read non-double (%s) for targetNum!\n", 
              tok);
            return -1;
        }
        thee->targetRes = tf;
        thee->settargetRes = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "maxsolve") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "parseFE:  read non-int (%s) for maxsolve!\n", tok);
            return -1;
        }
        thee->maxsolve = ti;
        thee->setmaxsolve = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "maxvert") == 0) {
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "parseFE:  read non-int (%s) for maxvert!\n", tok);
            return -1;
        }
        thee->maxvert = ti;
        thee->setmaxvert = 1;
        return 1;
    }


    return 0;

VERROR1:
    Vnm_print(2, "parseFE:  ran out of tokens!\n");
    return -1;

}
