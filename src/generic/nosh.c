/**
 *  @file    nosh.c
 *  @ingroup NOsh
 *  @author  Nathan Baker
 *  @brief   Class NOsh methods
 *  @version $Id$
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
 * Copyright (c) 2002.  Washington University in St. Louis.
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


#include "apbscfg.h"
#include "apbs/nosh.h"
#include "apbs/vstring.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Private method declaration
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int NOsh_parseREAD(NOsh *thee, Vio *sock);
VPRIVATE int NOsh_parsePRINT(NOsh *thee, Vio *sock);
VPRIVATE int NOsh_parseELEC(NOsh *thee, Vio *sock);
VPRIVATE int NOsh_parseFEM(NOsh *thee, Vio *sock, FEMparm *parm);
VEXTERNC int NOsh_parseMG(NOsh *thee, Vio *sock, int type);

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_NOSH)

    VPUBLIC char* NOsh_getMolpath(NOsh *thee, int imol) {
        VASSERT(thee != VNULL);
        VASSERT(imol < thee->nmol);
        return thee->molpath[imol];
    }
    VPUBLIC char* NOsh_getDielXpath(NOsh *thee, int imol) {
        VASSERT(thee != VNULL);
        VASSERT(imol < thee->nmol);
        return thee->dielXpath[imol];
    }
    VPUBLIC char* NOsh_getDielYpath(NOsh *thee, int imol) {
        VASSERT(thee != VNULL);
        VASSERT(imol < thee->nmol);
        return thee->dielYpath[imol];
    }
    VPUBLIC char* NOsh_getDielZpath(NOsh *thee, int imol) {
        VASSERT(thee != VNULL);
        VASSERT(imol < thee->nmol);
        return thee->dielZpath[imol];
    }
    VPUBLIC char* NOsh_getKappapath(NOsh *thee, int imol) {
        VASSERT(thee != VNULL);
        VASSERT(imol < thee->nmol);
        return thee->kappapath[imol];
    }
    VPUBLIC char* NOsh_getChargepath(NOsh *thee, int imol) {
        VASSERT(thee != VNULL);
        VASSERT(imol < thee->nmol);
        return thee->chargepath[imol];
    }
    VPUBLIC NOsh_calc* NOsh_getCalc(NOsh *thee, int icalc) {
        VASSERT(thee != VNULL);
        VASSERT(icalc < thee->ncalc);
        return &(thee->calc[icalc]);
    }
    VPUBLIC int NOsh_getDielfmt(NOsh *thee, int i) {
        VASSERT(thee != VNULL);
        VASSERT(i < thee->ndiel);
        return (thee->dielfmt[i]);
    }
    VPUBLIC int NOsh_getKappafmt(NOsh *thee, int i) {
        VASSERT(thee != VNULL);
        VASSERT(i < thee->nkappa);
        return (thee->kappafmt[i]);
    }
    VPUBLIC int NOsh_getChargefmt(NOsh *thee, int i) {
        VASSERT(thee != VNULL);
        VASSERT(i < thee->ncharge);
        return (thee->chargefmt[i]);
    }

#endif /* if !defined(VINLINE_NOSH) */

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

VPUBLIC int NOsh_printWhat(NOsh *thee, int iprint) {
    VASSERT(thee != VNULL);
    VASSERT(iprint < thee->nprint);
    return thee->printwhat[iprint];
}
VPUBLIC int NOsh_printNarg(NOsh *thee, int iprint) {
    VASSERT(thee != VNULL);
	VASSERT(iprint < thee->nprint);
    return thee->printnarg[iprint];
}
VPUBLIC int NOsh_elec2calc(NOsh *thee, int icalc) {
    VASSERT(thee != VNULL);
	VASSERT(icalc < thee->ncalc);
    return thee->elec2calc[icalc];
}
VPUBLIC int NOsh_printOp(NOsh *thee, int iprint, int iarg) {
    VASSERT(thee != VNULL);
    VASSERT(iprint < thee->nprint);
    VASSERT(iarg < thee->printnarg[iprint]);
    return thee->printop[iprint][iarg];
}
VPUBLIC int NOsh_printCalc(NOsh *thee, int iprint, int iarg) {
    VASSERT(thee != VNULL);
    VASSERT(iprint < thee->nprint);
    VASSERT(iarg < thee->printnarg[iprint]);
    return thee->printcalc[iprint][iarg];
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:   NOsh_ctor
//
// Author:    Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC NOsh* NOsh_ctor(int rank, int size) {

    /* Set up the structure */
    NOsh *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(NOsh) );
    VASSERT( thee != VNULL);
    VASSERT( NOsh_ctor2(thee, rank, size) );

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_ctor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int NOsh_ctor2(NOsh *thee, int rank, int size) {

    int i;

    if (thee == VNULL) return 0;

    thee->proc_rank = rank;
    thee->proc_size = size;
 
    thee->ispara = 0;
    thee->parsed = 0;
    thee->ncalc = 0;
    thee->nmol = 0;
    thee->ncharge = 0;
    thee->ndiel = 0;
    thee->nkappa = 0;
    thee->nprint = 0;

    for (i=0; i<NOSH_MAXCALC; i++) {
        thee->calc[i].mgparm = VNULL;
        thee->calc[i].femparm = VNULL;
        thee->calc[i].calctype = -1;
    }

    return 1; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_dtor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void NOsh_dtor(NOsh **thee) {
    if ((*thee) != VNULL) {
        NOsh_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(NOsh), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_dtor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void NOsh_dtor2(NOsh *thee) { 
   
    int i;

    if (thee != VNULL) {
        for (i=0; i<NOSH_MAXCALC; i++) {
            PBEparm_dtor(&(thee->calc[i].pbeparm));
            MGparm_dtor(&(thee->calc[i].mgparm));
            FEMparm_dtor(&(thee->calc[i].femparm));
        }
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_parseFile
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int NOsh_parseFile(NOsh *thee, char *filename) {

    Vio *sock;
    int rc;

    sock = Vio_ctor("FILE", "ASC", VNULL, filename, "r");
    rc = NOsh_parse(thee, sock);
    Vio_dtor(&sock);

    return rc;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_parse
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int NOsh_parse(NOsh *thee, Vio *sock) {
 
    char *MCwhiteChars = " =,;\t\n";
    char *MCcommChars  = "#%";
    char tok[VMAX_BUFSIZE];

    if (thee == VNULL) {
        Vnm_print(2, "NOsh_parse:  Got NULL thee!\n");
        return 0;
    }

    if (sock == VNULL) {
        Vnm_print(2, "NOsh_parse:  Got pointer to NULL socket!\n");
        return 0;
    } 

    if (thee->parsed) {
        Vnm_print(2, "NOsh_parse:  Already parsed an input file!\n");
        return 0;
    }

    if (Vio_accept(sock, 0) < 0) {
        Vnm_print(2, "NOsh_parse:  Problem reading from socket!\n");
        return 0;
    }

    /* Set up the whitespace and comment character definitions */
    Vio_setWhiteChars(sock, MCwhiteChars);
    Vio_setCommChars(sock, MCcommChars);

    /* We parse the file until we run out of tokens */
    Vnm_print(0, "NOsh_parse:  Starting file parsing...\n");
    while (Vio_scanf(sock, "%s", tok) == 1) {
        /* At the highest level, we look for keywords that indicate functions
         * like:
         *  read => Read in a molecule file
         *  elec => Do an electrostatics calculation
         * These cause the code to go to a lower-level parser routine which
         * handles keywords specific to the particular function.  Each
         * lower-level parser routine then returns when it hits the "end"
         * keyword.  Due to this simple layout, no nesting of these "function"
         * sections is allowed.
         */
        if (Vstring_strcasecmp(tok, "read") == 0) {
            Vnm_print(0, "NOsh: Parsing READ section\n");
            if (!NOsh_parseREAD(thee, sock)) return 0;
            Vnm_print(0, "NOsh: Done parsing READ section (nmol=%d, ndiel=%d, \
nkappa=%d, ncharge=%d)\n", thee->nmol, thee->ndiel, thee->nkappa, thee->ncharge);
        } else if (Vstring_strcasecmp(tok, "print") == 0) {
            Vnm_print(0, "NOsh: Parsing PRINT section\n");
            if (!NOsh_parsePRINT(thee, sock)) return 0;
            Vnm_print(0, "NOsh: Done parsing PRINT section\n");
        } else if (Vstring_strcasecmp(tok, "elec") == 0) {
            Vnm_print(0, "NOsh: Parsing ELEC section\n");
            if (!NOsh_parseELEC(thee, sock)) return 0;
            Vnm_print(0, "NOsh: Done parsing ELEC section (ncalc = %d)\n",
              thee->ncalc);
        } else if (Vstring_strcasecmp(tok, "quit") == 0) {
            Vnm_print(0, "NOsh: Done parsing file (got QUIT)\n");
            break;
        } else {
            Vnm_print(2, "NOsh_parse: Ignoring undefined keyword %s!\n", tok);
        }
    }

    thee->parsed = 1;
    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_parseREAD
//
// Purpose:  Parse an input file READ section
//
// Returns:  1 if successful, 0 otherwise
//
// Notes:    Should only be called from NOsh_parse()
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int NOsh_parseREAD(NOsh *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    int molfmt, dielfmt, chargefmt, kappafmt;

    if (thee == VNULL) {
        Vnm_print(2, "NOsh_parse:  Got NULL thee!\n");
        return 0;
    }

    if (sock == VNULL) {
        Vnm_print(2, "NOsh_parse:  Got pointer to NULL socket!\n");
        return 0;
    } 

    if (thee->parsed) {
        Vnm_print(2, "NOsh_parse:  Already parsed an input file!\n");
        return 0;
    }

    /* Read until we run out of tokens (bad) or hit the "END" keyword (good) */
    while (Vio_scanf(sock, "%s", tok) == 1) {

        if (Vstring_strcasecmp(tok, "end") == 0) {
            Vnm_print(0, "NOsh: Done parsing READ section\n");
            return 1;
        } else if (Vstring_strcasecmp(tok, "mol") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (Vstring_strcasecmp(tok, "pqr") == 0) {
                molfmt = 0;
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                Vnm_print(0, "NOsh: Storing molecule %d path %s\n", 
                  thee->nmol, tok);
                thee->molfmt[thee->nmol] = molfmt;
                strncpy(thee->molpath[thee->nmol], tok, VMAX_ARGLEN);
                (thee->nmol)++;
            } else {
                Vnm_print(2, "NOsh_parseREAD:  Ignoring undefined format \
%s!\n", tok);
            } 
        } else if (Vstring_strcasecmp(tok, "diel") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (Vstring_strcasecmp(tok, "dx") == 0) {
                dielfmt = 0;
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                Vnm_print(0, "NOsh: Storing x-shifted dielectric map %d path \
%s\n", thee->ndiel, tok);
                strncpy(thee->dielXpath[thee->ndiel], tok, VMAX_ARGLEN);
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                Vnm_print(0, "NOsh: Storing y-shifted dielectric map %d path \
%s\n", thee->ndiel, tok);
                strncpy(thee->dielYpath[thee->ndiel], tok, VMAX_ARGLEN);
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                Vnm_print(0, "NOsh: Storing z-shifted dielectric map %d path \
%s\n", thee->ndiel, tok);
                strncpy(thee->dielZpath[thee->ndiel], tok, VMAX_ARGLEN);
                thee->dielfmt[thee->ndiel] = dielfmt;
                (thee->ndiel)++;
            } else { 
                Vnm_print(2, "NOsh_parseREAD:  Ignoring undefined format \
%s!\n", tok);
            } 
        } else if (Vstring_strcasecmp(tok, "kappa") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (Vstring_strcasecmp(tok, "dx") == 0) {
                kappafmt = 0;
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                Vnm_print(0, "NOsh: Storing kappa map %d path %s\n",
                  thee->nkappa, tok);
                thee->kappafmt[thee->nkappa] = kappafmt;
                strncpy(thee->kappapath[thee->nkappa], tok, VMAX_ARGLEN);
                (thee->nkappa)++;
            } else {
                Vnm_print(2, "NOsh_parseREAD:  Ignoring undefined format \
%s!\n", tok);
            }
        } else if (Vstring_strcasecmp(tok, "charge") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (Vstring_strcasecmp(tok, "dx") == 0) {
                chargefmt = 0;
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                Vnm_print(0, "NOsh: Storing charge map %d path %s\n",
                  thee->ncharge, tok);
                thee->chargefmt[thee->ncharge] = chargefmt;
                strncpy(thee->chargepath[thee->ncharge], tok, VMAX_ARGLEN);
                (thee->ncharge)++;
            } else {
                Vnm_print(2, "NOsh_parseREAD:  Ignoring undefined format \
%s!\n", tok);
            }
        } else {
            Vnm_print(2, "NOsh_parseREAD:  Ignoring undefined keyword %s!\n", 
              tok);
        }
    }

    /* We ran out of tokens! */
    VJMPERR1(0);

    VERROR1:
        Vnm_print(2, "NOsh_parseREAD:  Ran out of tokens while parsing READ \
section!\n");
        return 0;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_parsePRINT
//
// Purpose:  Parse an input file PRINT section
//
// Returns:  1 if successful, 0 otherwise
//
// Notes:    Should only be called from NOsh_parse()
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int NOsh_parsePRINT(NOsh *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    int ti, idx, expect;

    if (thee == VNULL) {
        Vnm_print(2, "NOsh_parse:  Got NULL thee!\n");
        return 0;
    }

    if (sock == VNULL) {
        Vnm_print(2, "NOsh_parse:  Got pointer to NULL socket!\n");
        return 0;
    } 

    if (thee->parsed) {
        Vnm_print(2, "NOsh_parse:  Already parsed an input file!\n");
        return 0;
    }

    idx = thee->nprint;
    if (thee->nprint >= NOSH_MAXPRINT) {
        Vnm_print(2, "NOsh_parse:  Exceeded max number (%d) of PRINT \
sections\n", 
          NOSH_MAXPRINT);
        return 0;
    }

   
    /* The first thing we read is the thing we want to print */ 
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (Vstring_strcasecmp(tok, "energy") == 0) {
        thee->printwhat[idx] = 0;
        thee->printnarg[idx] = 0;
    } else if (Vstring_strcasecmp(tok, "force") == 0) {
        thee->printwhat[idx] = 1;
        thee->printnarg[idx] = 0;
    } else {
        Vnm_print(2, "NOsh_parsePRINT:  Undefined keyword %s while parsing \
PRINT section!\n",
          tok);
        return 0;
    }

    expect = 0;   /* We first expect a calculation ID (0) then an op (1) */

    /* Read until we run out of tokens (bad) or hit the "END" keyword (good) */
    while (Vio_scanf(sock, "%s", tok) == 1) {

        /* The next thing we read is either END or an ARG OP ARG statement */
        if (Vstring_strcasecmp(tok, "end") == 0) {
            if (expect != 0) {
                (thee->nprint)++;
                (thee->printnarg[idx])++;
                Vnm_print(0, "NOsh: Done parsing PRINT section\n");
                return 1;
            } else {
                Vnm_print(2, "NOsh_parsePRINT:  Got premature END to PRINT!\n");
                return 0;
            }
        } else {
            /* Grab a calculation ID */
            if (sscanf(tok, "%d", &ti) == 1) {
                if (expect == 0) {
                    thee->printcalc[idx][thee->printnarg[idx]] = ti;
                    expect = 1;
                } else {
                    Vnm_print(2, "NOsh_parsePRINT:  Syntax error in PRINT \
section while reading %s!\n", tok);
                    return 0;
                }
            /* Grab addition operation */
            } else if (Vstring_strcasecmp(tok, "+") == 0) {
                if (expect == 1) {
                    thee->printop[idx][thee->printnarg[idx]] = 0;
                    (thee->printnarg[idx])++;
                    expect = 0;
                    if (thee->printnarg[idx] >= NOSH_MAXPOP) {
                        Vnm_print(2, "NOsh_parsePRINT:  Exceeded max number \
(%d) of arguments for PRINT section!\n", 
                          NOSH_MAXPOP);
                        return 0;
                    }
                } else {
                    Vnm_print(2, "NOsh_parsePRINT:  Syntax error in PRINT \
section while reading %s!\n", tok);
                    return 0;
                }
            /* Grab subtraction operation */
            } else if (Vstring_strcasecmp(tok, "-") == 0) {
                if (expect == 1) {
                    thee->printop[idx][thee->printnarg[idx]] = 1;
                    (thee->printnarg[idx])++;
                    expect = 0;
                    if (thee->printnarg[idx] >= NOSH_MAXPOP) {
                        Vnm_print(2, "NOsh_parseREAD:  Exceeded max number \
(%d) of arguments for PRINT section!\n", 
                          NOSH_MAXPOP);
                        return 0;
                    }
                } else {
                    Vnm_print(2, "NOsh_parsePRINT:  Syntax error in PRINT \
section while reading %s!\n", tok);
                    return 0;
                }
            /* Got bad operation */
            } else {
                Vnm_print(2, "NOsh_parsePRINT:  Undefined keyword %s while \
parsing PRINT section!\n",
                  tok);
                return 0;
            } 
        } /* end parse token */

    } /* end while */

    VJMPERR1(0);

    /* We ran out of tokens! */
    VERROR1:
       Vnm_print(2, "NOsh_parsePRINT:  Ran out of tokens while parsing PRINT \
section!\n");
       return 0;

}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_parseELEC
//
// Purpose:  Parse an input file ELEC section
//
// Returns:  1 if successful, 0 otherwise
//
// Notes:    Should only be called from NOsh_parse()
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int NOsh_parseELEC(NOsh *thee, Vio *sock) {
 
    MGparm *tmgparms[NOSH_MAXCALC];
    int i;

    char tok[VMAX_BUFSIZE];

    if (thee == VNULL) {
        Vnm_print(2, "NOsh_parseELEC:  Got NULL thee!\n");
        return 0;
    }

    if (sock == VNULL) {
        Vnm_print(2, "NOsh_parseELEC:  Got pointer to NULL socket!\n");
        return 0;
    } 

    if (thee->parsed) {
        Vnm_print(2, "NOsh_parseELEC:  Already parsed an input file!\n");
        return 0;
    }

    for (i=0; i<NOSH_MAXCALC; i++) tmgparms[i] = VNULL;

    /* Update the ELEC statement number */
    if (thee->ncalc >= NOSH_MAXCALC) {
        Vnm_print(2, "NOsh:  Too many electrostatics calculations in this \
run!\n");
        Vnm_print(2, "NOsh:  Current max is %d; ignoring this calculation\n", 
          NOSH_MAXCALC);
        return 1;
    } else (thee->nelec)++;

    /* The next token HAS to be the method */
    if (Vio_scanf(sock, "%s", tok) == 1) {
        if (Vstring_strcasecmp(tok, "mg-manual") == 0) {
            return NOsh_parseMG(thee, sock, 0);
        } else if (Vstring_strcasecmp(tok, "mg-auto") == 0) {
            return NOsh_parseMG(thee, sock, 1);
        } else if (Vstring_strcasecmp(tok, "mg-para") == 0) {
            return NOsh_parseMG(thee, sock, 2);
        } else if (Vstring_strcasecmp(tok, "mg-dummy") == 0) {
            return NOsh_parseMG(thee, sock, 3);
        } else if (Vstring_strcasecmp(tok, "fem") == 0) {
            /* Check to see if he have any room left for this type of
             * calculation, if so: set the calculation type, update the number
             * of calculations of this type, and parse the rest of the section
             */
            if (thee->ncalc >= NOSH_MAXCALC) {
                Vnm_print(2, "NOsh:  Too many calculations in this run!\n");
                Vnm_print(2, "NOsh:  Current max is %d; ignoring this calculation\n",
                  NOSH_MAXCALC);
                return 1;
            }
            (thee->ncalc)++;
            thee->calc[thee->ncalc - 1].calctype = 1;
            Vnm_print(0, "NOsh: Parsing parameters for FEM calculation #%d\n",
              thee->ncalc);
            thee->calc[thee->ncalc-1].femparm = FEMparm_ctor();
            thee->elec2calc[thee->nelec-1] = thee->ncalc-1;
            return NOsh_parseFEM(thee,sock,thee->calc[thee->ncalc-1].femparm);
        } else {
            Vnm_print(2, "NOsh_parseELEC: The method (\"mg\" or \"fem\") must be the first keyword in the ELEC section\n");
            return 0;
        }
    } else {
        Vnm_print(2, "NOsh_parseELEC:  Ran out of tokens while reading ELEC section!\n");
        return 0;
    } 

    /* We ran out of tokens! */
    Vnm_print(2, "NOsh_parseELEC:  Ran out of tokens while parsing ELEC section!\n");
    return 0;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_parseFEM
//
// Purpose:  Parse an input file ELEC section for the FEM method
//
// Returns:  1 if successful, 0 otherwise
//
// Notes:    Should only be called from NOsh_parse()
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int NOsh_parseFEM(NOsh *thee, Vio *sock, FEMparm *parm) {

    if (thee == VNULL) {
        Vnm_print(2, "NOsh_parseFEM:  Got NULL thee!\n");
        return 0;
    }

    if (parm == VNULL) {
        Vnm_print(2, "NOsh_parseFEM:  Got NULL parm!\n");
        return 0;
    }

    if (sock == VNULL) {
        Vnm_print(2, "NOsh_parseFEM:  Got pointer to NULL socket!\n");
        return 0;
    }

    if (thee->parsed) {
        Vnm_print(2, "NOsh_parseFEM:  Already parsed an input file!\n");
        return 0;
    }

    Vnm_print(2, "NOsh_parseFEM:  FEM not availble yet; igoring this section!\n");
    return 1;
}
