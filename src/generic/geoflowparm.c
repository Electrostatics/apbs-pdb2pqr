/**
 *  @file    geoflowparm.c
 *  @ingroup GEOFLOWparm
 *  @author  Andrew Stevens
 *  @brief   Class GEOFLOWparm methods
 *  @version $Id$
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 *  Nathan A. Baker (nathan.baker@pnnl.gov)
 *  Pacific Northwest National Laboratory
 *
 *  Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2012 Battelle Memorial Institute. Developed at the
 * Pacific Northwest National Laboratory, operated by Battelle Memorial
 * Institute, Pacific Northwest Division for the U.S. Department of Energy.
 *
 * Portions Copyright (c) 2002-2010, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2010, Nathan A. Baker.
 * Portions Copyright (c) 1999-2002, The Regents of the University of
 * California.
 * Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * Neither the name of the developer nor the names of its contributors may be
 * used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#include "geoflowparm.h"

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_MGPARM)

#endif /* if !defined(VINLINE_MGPARM) */


VPUBLIC GEOFLOWparm* GEOFLOWparm_ctor(GEOFLOWparm_CalcType type) {

    /* Set up the structure */
    GEOFLOWparm *thee = VNULL;
    thee = (GEOFLOWparm*)Vmem_malloc(VNULL, 1, sizeof(GEOFLOWparm));
    VASSERT( thee != VNULL);
    VASSERT( GEOFLOWparm_ctor2(thee, type) == VRC_SUCCESS );

    return thee;
}

VPUBLIC Vrc_Codes GEOFLOWparm_ctor2(GEOFLOWparm *thee, GEOFLOWparm_CalcType type) {

    int i;

    if (thee == VNULL) return VRC_FAILURE;

    thee->parsed = 0;
    thee->type = type;

    /* *** GENERIC PARAMETERS *** */

    /* *** TYPE 0 PARAMETERS *** */
    thee->tree_order = 1;
    thee->settree_order = 0;
    thee->tree_n0 = 500;
    thee->settree_n0 = 0;
    thee->mac = 0.8;
    thee->setmac = 0;

    /* *** TYPE 1 & 2 PARAMETERS *** */

    /* *** TYPE 2 PARAMETERS *** */
    thee->nonlintype = 0;
    thee->setnonlintype = 0;

    /* *** Default parameters for TINKER *** */
    thee->chgs = VCM_CHARGE;

    return VRC_SUCCESS;
}

VPUBLIC void GEOFLOWparm_dtor(GEOFLOWparm **thee) {
    if ((*thee) != VNULL) {
        GEOFLOWparm_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(GEOFLOWparm), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void GEOFLOWparm_dtor2(GEOFLOWparm *thee) { ; }

VPUBLIC Vrc_Codes GEOFLOWparm_check(GEOFLOWparm *thee) {

    Vrc_Codes rc;
    int i, tdime[3], ti, tnlev[3], nlev;

    rc = VRC_SUCCESS;

    Vnm_print(0, "GEOFLOWparm_check:  checking GEOFLOWparm object of type %d.\n",
              thee->type);

    /* Check to see if we were even filled... */
    if (!thee->parsed) {
        Vnm_print(2, "GEOFLOWparm_check:  not filled!\n");
        return VRC_FAILURE;
    }


    /* Check type settings */
    if ((thee->type != GFCT_MANUAL) && (thee->type != GFCT_NONE)) {
         Vnm_print(2,"GEOFLOWparm_check: type not set");
         rc = VRC_FAILURE;
    }

    /* Check treecode setting*/
    if (thee->tree_order<1) {
        Vnm_print(2,"GEOFLOWparm_check: treecode order is less than 1");
        rc = VRC_FAILURE;
    }
    if (thee->tree_n0<1) {
        Vnm_print(2,"GEOFLOWparm_check: treecode leaf size is less than 1");
        rc = VRC_FAILURE;
    }
    if (thee->mac>1 || thee->mac <=0 ) {
        Vnm_print(2,"GEOFLOWparm_check: MAC criterion fails");
        rc = VRC_FAILURE;
    }

    return rc;
}

VPUBLIC void GEOFLOWparm_copy(GEOFLOWparm *thee, GEOFLOWparm *parm) {

    int i;

    VASSERT(thee != VNULL);
    VASSERT(parm != VNULL);


    thee->type = parm->type;
    thee->parsed = parm->parsed;

    /* *** GENERIC PARAMETERS *** */

    /* *** TYPE 0 PARMS *** */

    thee->tree_order = parm->tree_order;
    thee->settree_order = parm->settree_order;
    thee->tree_n0 = parm->tree_n0;
    thee->settree_n0 = parm->settree_n0;
    thee->mac = parm->mac;
    thee->setmac = parm->setmac;

    /* *** TYPE 1 & 2 PARMS *** */

    /* *** TYPE 2 PARMS *** */
    thee->nonlintype = parm->nonlintype;
    thee->setnonlintype = parm->setnonlintype;
}


VPRIVATE Vrc_Codes GEOFLOWparm_parseTREE_ORDER(GEOFLOWparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    int ti;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%d", &ti) == 0) {
        Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing TREE_ORDER \
keyword!\n", tok);
        return VRC_WARNING;
    } else if (ti <= 0) {
        Vnm_print(2, "parseGEOFLOW:  tree_order must be greater than 0!\n");
        return VRC_WARNING;
    } else thee->tree_order = ti;
    thee->settree_order = 1;
    return VRC_SUCCESS;

    VERROR1:
        Vnm_print(2, "parseGEOFLOW:  ran out of tokens!\n");
        return VRC_WARNING;
}


VPRIVATE Vrc_Codes GEOFLOWparm_parseTREE_N0(GEOFLOWparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    int ti;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%d", &ti) == 0) {
        Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing TREE_ORDER \
keyword!\n", tok);
        return VRC_WARNING;
    } else if (ti <= 0) {
        Vnm_print(2, "parseGEOFLOW:  tree_order must be greater than 0!\n");
        return VRC_WARNING;
    } else thee->tree_order = ti;
    thee->settree_order = 1;
    return VRC_SUCCESS;

    VERROR1:
        Vnm_print(2, "parseGEOFLOW:  ran out of tokens!\n");
        return VRC_WARNING;
}
VPRIVATE Vrc_Codes GEOFLOWparm_parseMAC(GEOFLOWparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing mac \
keyword!\n", tok);
        return VRC_WARNING;
    } else if (tf <= 0.0 || tf > 1.0) {
        Vnm_print(2, "parseGEOFLOW:  mac must be between 0 and 1!\n");
        return VRC_WARNING;
    } else thee->mac = tf;
    thee->setmac = 1;
    return VRC_SUCCESS;

    VERROR1:
        Vnm_print(2, "parseGEOFLOW:  ran out of tokens!\n");
        return VRC_WARNING;
}


VPUBLIC Vrc_Codes GEOFLOWparm_parseToken(GEOFLOWparm *thee, char tok[VMAX_BUFSIZE],
  Vio *sock) {

    if (thee == VNULL) {
        Vnm_print(2, "parseGEOFLOW:  got NULL thee!\n");
        return VRC_WARNING;
    }
    if (sock == VNULL) {
        Vnm_print(2, "parseGEOFLOW:  got NULL socket!\n");
        return VRC_WARNING;
    }

    Vnm_print(0, "GEOFLOWparm_parseToken:  trying %s...\n", tok);


    if (Vstring_strcasecmp(tok, "tree_order") == 0) {
        return GEOFLOWparm_parseTREE_ORDER(thee, sock);
    } else if (Vstring_strcasecmp(tok, "tree_n0") == 0) {
        return GEOFLOWparm_parseTREE_N0(thee, sock);
    } else if (Vstring_strcasecmp(tok, "mac") == 0) {
        return GEOFLOWparm_parseMAC(thee, sock);
    } else {
        Vnm_print(2, "parseGEOFLOW:  Unrecognized keyword (%s)!\n", tok);
        return VRC_WARNING;
    }

    return VRC_FAILURE;

}
