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
// File:     nosh.c
//
// Purpose:  Class NOsh: methods. 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "apbs/nosh.h"
#include "mgautoparm.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Private method declaration
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int NOsh_parseREAD(NOsh *thee, Vio *sock);
VPRIVATE int NOsh_parsePRINT(NOsh *thee, Vio *sock);
VPRIVATE int NOsh_parseELEC(NOsh *thee, Vio *sock);
VPRIVATE int NOsh_parseFEM(NOsh *thee, Vio *sock, FEMparm *parm);
VPRIVATE int NOsh_parseMGMANUAL(NOsh *thee, Vio *sock, MGparm *parm);
VPRIVATE int NOsh_parseMGAUTO(NOsh *thee, Vio *sock, int *nparm, 
  MGparm *parms[NOSH_MAXCALC]);

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_NOSH)

#endif /* if !defined(VINLINE_NOSH) */

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_ctor
//
// Author: Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC NOsh* NOsh_ctor() {

    /* Set up the structure */
    NOsh *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(NOsh) );
    VASSERT( thee != VNULL);
    VASSERT( NOsh_ctor2(thee) );

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_ctor2
//
// Purpose:  Construct the NOsh object
//
// Notes:    Constructor broken into two parts for FORTRAN users.
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int NOsh_ctor2(NOsh *thee) {

    int i;

    if (thee == VNULL) return 0;

    thee->parsed = 0;
    thee->ncalc = 0;
    thee->nmol = 0;
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
// Notes:    Since the grid manager and atom list were allocated outside of
//           the NOsh routines, they are not destroyed.
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
            if (thee->calc[i].mgparm != VNULL) 
              MGparm_dtor(&(thee->calc[i].mgparm));
            if (thee->calc[i].femparm != VNULL) 
              FEMparm_dtor(&(thee->calc[i].femparm));
        }
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_parse
//
// Purpose:  Parse an input file
//
// Returns:  1 if successful, 0 otherwise
//
// Notes:    Currently uses strcasecmp(), which I think is not ANSI C compliant
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
        if (strcasecmp(tok, "read") == 0) {
            Vnm_print(0, "NOsh: Parsing READ section\n");
            if (!NOsh_parseREAD(thee, sock)) return 0;
            Vnm_print(0, "NOsh: Done parsing READ section (nmol = %d)\n",
              thee->nmol);
        } else if (strcasecmp(tok, "print") == 0) {
            Vnm_print(0, "NOsh: Parsing PRINT section\n");
            if (!NOsh_parsePRINT(thee, sock)) return 0;
            Vnm_print(0, "NOsh: Done parsing PRINT section\n");
        } else if (strcasecmp(tok, "elec") == 0) {
            Vnm_print(0, "NOsh: Parsing ELEC section\n");
            if (!NOsh_parseELEC(thee, sock)) return 0;
            Vnm_print(0, "NOsh: Done parsing ELEC section (ncalc = %d)\n",
              thee->ncalc);
        } else if (strcasecmp(tok, "quit") == 0) {
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
// Notes:    Currently uses strcasecmp(), which I think is not ANSI C compliant
//           Should only be called from NOsh_parse()
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int NOsh_parseREAD(NOsh *thee, Vio *sock) {

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

    /* Read until we run out of tokens (bad) or hit the "END" keyword (good) */
    while (Vio_scanf(sock, "%s", tok) == 1) {

        if (strcasecmp(tok, "end") == 0) {
            Vnm_print(0, "NOsh: Done parsing READ section\n");
            return 1;
        } else if (strcasecmp(tok, "file") == 0) {
            if (Vio_scanf(sock, "%s", tok) == 1) {
                Vnm_print(0, "NOsh: Storing molecule %d path %s\n", 
                  thee->nmol, tok);
                strncpy(thee->molpath[thee->nmol], tok, VMAX_ARGLEN);
                (thee->nmol)++;
            } else break;
        } else {
            Vnm_print(2, "NOsh_parseREAD:  Ignoring undefined keyword %s!\n", 
              tok);
        }

    }

    /* We ran out of tokens! */
    Vnm_print(2, "NOsh_parseREAD:  Ran out of tokens while parsing READ section!\n");
    return 0;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_parsePRINT
//
// Purpose:  Parse an input file PRINT section
//
// Returns:  1 if successful, 0 otherwise
//
// Notes:    Currently uses strcasecmp(), which I think is not ANSI C compliant
//           Should only be called from NOsh_parse()
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
    if (strcasecmp(tok, "energy") == 0) {
        thee->printwhat[idx] = 0;
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
        if (strcasecmp(tok, "end") == 0) {
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
            } else if (strcasecmp(tok, "+") == 0) {
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
            } else if (strcasecmp(tok, "-") == 0) {
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
// Notes:    Currently uses strcasecmp(), which I think is not ANSI C compliant
//           Should only be called from NOsh_parse()
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int NOsh_parseELEC(NOsh *thee, Vio *sock) {
 
    MGparm *tmgparms[NOSH_MAXCALC];
    int i, nparm;

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
        if ((strcasecmp(tok, "mg") == 0) || 
            (strcasecmp(tok, "mg-manual") == 0)) {
            if (strcasecmp(tok, "mg") == 0) { 
                Vnm_print(2, "NOsh:  The MG keyword is deprecated.  Please use\n");
                Vnm_print(2, "NOsh:  either MG-MANUAL or MG-AUTO.  I'm assuming\n");
                Vnm_print(2, "NOsh:  you meant MG-MANUAL here.\n");
            }
	    /* Check to see if he have any room left for this type of
             * calculation, if so: set the calculation type, update the number
             * of calculations of this type, and parse the rest of the section
             */
            if (thee->ncalc >= NOSH_MAXCALC) {
                Vnm_print(2, "NOsh:  Too many calculations in this run!\n");
                Vnm_print(2, "NOsh:  Current max is %d; ignoring this \
calculation\n",
                  NOSH_MAXCALC);
                return 1;
            }
            (thee->ncalc)++;
            thee->calc[thee->ncalc-1].calctype = 0;
            Vnm_print(0, "NOsh: Parsing parameters for MG calculation #%d\n",
              thee->ncalc);
            thee->calc[thee->ncalc-1].mgparm = MGparm_ctor();
            thee->elec2calc[thee->nelec-1] = thee->ncalc-1;
            return NOsh_parseMGMANUAL(thee, sock, 
              thee->calc[thee->ncalc-1].mgparm);
        } else if (strcasecmp(tok, "mg-auto") == 0) {
            Vnm_print(0, "NOsh: Parsing parameters for MG-AUTO calculation.\n");
            if (!NOsh_parseMGAUTO(thee, sock, &nparm, tmgparms)) return 0;
            if ((thee->ncalc + nparm) >= NOSH_MAXCALC) {
                Vnm_print(2, "NOsh:  Foucsing requires too many multigrid \
electrostatics calculations in this run!\n");
                Vnm_print(2, "NOsh:  Current max is %d; ignoring this \
calculation\n",
                  NOSH_MAXCALC);
                return 1;
            }
            for (i=0; i<nparm; i++) {
                thee->calc[thee->ncalc+i].mgparm = tmgparms[i];
                thee->calc[thee->ncalc+i].calctype = 0;
            }
            /* Setup the map from ELEC statement to actual calculation */
            thee->elec2calc[thee->nelec] = thee->ncalc + nparm - 1;
            thee->ncalc += nparm;
        } else if (strcasecmp(tok, "fem") == 0) {
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
// Notes:    Currently uses strcasecmp(), which I think is not ANSI C compliant
//           Should only be called from NOsh_parse()
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int NOsh_parseFEM(NOsh *thee, Vio *sock, FEMparm *parm) {

    char tok[VMAX_BUFSIZE];

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

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_parseMGMANUAL
//
// Purpose:  Parse an input file ELEC section for the MG method
//
// Returns:  1 if successful, 0 otherwise
//
// Notes:    Currently uses strcasecmp(), which I think is not ANSI C compliant
//           Should only be called from NOsh_parse()
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int NOsh_parseMGMANUAL(NOsh *thee, Vio *sock, MGparm *parm) {

    char tok[VMAX_BUFSIZE];
    double tf;
    int ti;

    if (thee == VNULL) {
        Vnm_print(2, "NOsh:  Got NULL thee!\n");
        return 0;
    }

    if (parm == VNULL) {
        Vnm_print(2, "NOsh:  Got NULL parm!\n");
        return 0;
    }

    if (sock == VNULL) {
        Vnm_print(2, "NOsh:  Got pointer to NULL socket!\n");
        return 0;
    }

    if (thee->parsed) {
        Vnm_print(2, "NOsh:  Already parsed an input file!\n");
        return 0;
    }

    /* Here we go... */
    while (Vio_scanf(sock, "%s", tok) == 1) {
        if (strcasecmp(tok, "end") == 0) {
            /* Check to see that everything was set */
            parm->parsed = 1;
            if (!MGparm_check(parm)) { 
                Vnm_print(2, "NOsh:  MG parameters not set correctly!\n");
                VJMPERR2(0);
            }
            Vnm_print(0, "NOsh:  Done parsing ELEC section\n");
            return 1;
        /* Read grid dimensions */
        } else if (strcasecmp(tok, "dime") == 0) {
            /* Read the number of grid points */
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing DIME keyword!\n",
                  tok);
                VJMPERR2(0);
            } else parm->dime[0] = ti;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing DIME keyword!\n",
                  tok);
                VJMPERR2(0);
            } else parm->dime[1] = ti;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing DIME keyword!\n",
                  tok);
                VJMPERR2(0);
            } else parm->dime[2] = ti;
            parm->setdime = 1;
        /* Read number of levels in hierarchy */
        } else if (strcasecmp(tok, "nlev") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing NLEV keyword!\n",
                  tok);
                VJMPERR2(0);
            } else parm->nlev = ti;
            parm->setnlev = 1;
        } else if (strcasecmp(tok, "grid") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GRID keyword!\n",
                  tok);
                VJMPERR2(0);
            } else parm->grid[0] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GRID keyword!\n",
                  tok);
                VJMPERR2(0);
            } else parm->grid[1] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GRID keyword!\n",
                  tok);
                VJMPERR2(0);
            } else parm->grid[2] = tf;
            parm->setgrid = 1;
        /* Grid length keyword */
        } else if (strcasecmp(tok, "glen") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN keyword!\n",
                  tok);
                VJMPERR2(0);
            } else parm->glen[0] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN keyword!\n",
                  tok);
                VJMPERR2(0);
            } else parm->glen[1] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN keyword!\n",
                  tok);
                VJMPERR2(0);
            } else parm->glen[2] = tf;
            parm->setglen = 1;
        /* Grid center keyword */
        } else if (strcasecmp(tok, "gcent") == 0) {
            /* If the next token isn't a float, it probably means we want to
             * center on a molecule */
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                if (strcasecmp(tok, "mol") == 0) {
                    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                    if (sscanf(tok, "%d", &ti) == 0) {
                        Vnm_print(2, "NOsh:  Read non-int (%s) while parsing GCENT MOL keyword!\n",
                          tok);
                        VJMPERR2(0);
                    } else {
                        parm->cmeth = 1;
                        parm->centmol = ti;
                    }
                } else {
                    Vnm_print(2, "NOsh:  Unexpected keyword (%s) while parsing GCENT!\n",
                      tok);
                    VJMPERR2(0);
                }
            } else {
                parm->center[0] = tf;
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%lf", &tf) == 0) {
                    Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GCENT keyword!\n",
                      tok);
                    VJMPERR2(0);
                } 
                parm->center[1] = tf;
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%lf", &tf) == 0) {
                    Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GCENT keyword!\n",
                      tok);
                    VJMPERR2(0);
                } 
                parm->center[2] = tf;
            }   
            parm->setgcent = 1;
        /* Read mol ID */
        } else if (strcasecmp(tok, "mol") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-int (%s) while parsing MOL keyword!\n",
                  tok);
                VJMPERR2(0);
            } 
            parm->molid = ti;
            parm->setmolid = 1;
        /* Linearized vs. nonlinear PBE */
        } else if (strcasecmp(tok, "lpbe") == 0) {
            parm->nonlin = 0;
            parm->setnonlin = 1;
        } else if (strcasecmp(tok, "npbe") == 0) {
            parm->nonlin = 1;
            parm->setnonlin = 1;
        /* Boundary condition flag */
        } else if (strcasecmp(tok, "bcfl") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-int (%s) while parsing BCFL keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->bcfl = ti;
            parm->setbcfl = 1;
        /* Ions */
        } else if (strcasecmp(tok, "ion") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->ionq[parm->nion] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->ionc[parm->nion] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->ionr[parm->nion] = tf;
            parm->setion[parm->nion] = 1;
            (parm->nion)++;
            parm->setnion = 1;
        /* Solute dielectric */
        } else if (strcasecmp(tok, "pdie") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing PDIE keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->pdie = tf;
            parm->setpdie = 1;
        /* Solvent dielectric */
        } else if (strcasecmp(tok, "sdie") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SDIE keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->sdie = tf;
            parm->setsdie = 1;
        /* Surface definition methods */
        } else if (strcasecmp(tok, "srfm") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-int (%s) while parsing SRFM keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->srfm = ti;
            parm->setsrfm = 1;
        /* Solvent radius */
        } else if (strcasecmp(tok, "srad") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SRAD keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->srad = tf;
            parm->setsrad = 1;
        /* Spline window */
        } else if (strcasecmp(tok, "swin") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SWIN keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->swin = tf;
            parm->setswin = 1;
        /* Temperature */
        } else if (strcasecmp(tok, "temp") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing TEMP keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->temp = tf;
            parm->settemp = 1;
        /* Surface tension (apolar) */
        } else if (strcasecmp(tok, "gamma") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GAMMA keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->gamma = tf;
            parm->setgamma = 1;
        /* Energy writing */
        } else if (strcasecmp(tok, "calcenergy") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing WRITEENERGY keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->calcenergy = ti;
            parm->setcalcenergy = 1;
        /* Force writing */
        } else if (strcasecmp(tok, "calcforce") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing WRITEFORCE keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->calcforce = ti;
            parm->setcalcforce = 1;
        /* Potential writing */
        } else if (strcasecmp(tok, "writepot") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing WRITEPOT keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->writepot = ti;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (strcasecmp(tok, "dx") == 0) {
                parm->writepotfmt = 0;
            } else if (strcasecmp(tok, "avs") == 0) {
                parm->writepotfmt = 1;
            } else if (strcasecmp(tok, "uhbd") == 0) {
                parm->writepotfmt = 2;
            } else {
                Vnm_print(2, "NOsh:  Invalid format (%s) while parsing WRITEPOT keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            strncpy(parm->writepotstem, tok, VMAX_ARGLEN); 
            parm->setwritepot = 1;
         /* Accessibility writing */
        } else if (strcasecmp(tok, "writeacc") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing WRITEACC keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            parm->writeacc = ti;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (strcasecmp(tok, "dx") == 0) {
                parm->writeaccfmt = 0;
            } else if (strcasecmp(tok, "avs") == 0) {
                parm->writeaccfmt = 1;
            } else if (strcasecmp(tok, "uhbd") == 0) {
                parm->writeaccfmt = 2;
            } else {
                Vnm_print(2, "NOsh:  Invalid format (%s) while parsing \
WRITEPOT keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            strncpy(parm->writeaccstem, tok, VMAX_ARGLEN);
            parm->setwriteacc = 1;
        } else {
            Vnm_print(2, "NOsh:  Ignoring unknown keyword (%s) while parsing ELEC section!\n",
              tok);
        }
    }
  
    /* Ran out of tokens? */
    Vnm_print(2, "NOsh:  Ran out of tokens while parsing ELEC section!\n");
    VJMPERR2(0);

    VERROR1:
       Vnm_print(2, "NOsh:  Ran out of tokens while parsing ELEC \
section!\n");
       VJMPERR2(0);

    /* Default garbage collection handler during an error */
    VERROR2:
        return 0;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_parseMGAUTO
//
// Purpose:  Parse an input file ELEC section for the MG method
//
// Returns:  1 if successful, 0 otherwise
//
// Args:     sock    Socket (file) object
//           nparm   Number of MGparm objects created
//           parms   Array (length nparm) of pointers to MGparm objects
//
// Notes:    Currently uses strcasecmp(), which I think is not ANSI C compliant
//           Should only be called from NOsh_parse()
//           The user is responsible for destroying the MGparm objects when
//           done
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int NOsh_parseMGAUTO(NOsh *thee, Vio *sock, int *nparm, 
  MGparm *parms[NOSH_MAXCALC]) {

    MGparm *cparm, *fparm;
    MGAUTOparm *mgauto;
    char tok[VMAX_BUFSIZE];
    double tf;
    int i, ti;

    if (thee == VNULL) {
        Vnm_print(2, "NOsh:  Got NULL thee!\n");
        return 0;
    }

    if (parms != VNULL) {
        Vnm_print(2, "NOsh:  Destroying existing parms!\n");
        return 0;
    }

    if (sock == VNULL) {
        Vnm_print(2, "NOsh:  Got pointer to NULL socket!\n");
        return 0;
    }

    if (thee->parsed) {
        Vnm_print(2, "NOsh:  Already parsed an input file!\n");
        return 0;
    }

    /* Construct the coarse and fine parameter objects */
    cparm = MGparm_ctor();
    fparm = MGparm_ctor();

    /* Here we go... */
    while (Vio_scanf(sock, "%s", tok) == 1) {
        if (strcasecmp(tok, "end") == 0) {
            /* Check to see that everything was set */
            cparm->parsed = 1;
            fparm->parsed = 1;
            if (!MGparm_check(cparm) || !MGparm_check(fparm)) { 
                Vnm_print(2, "NOsh: MG parameters not set correctly!\n");
                VJMPERR2(0);
            }
            /* Build the real parameter objects */
            mgauto = MGAUTOparm_ctor(cparm, fparm);
            MGAUTOparm_build(mgauto, nparm, parms);
            Vnm_print(0, "NOsh:  Built %d parameter objects for focusing\n",
              nparm);
            MGAUTOparm_dtor(&mgauto);
            /* Destroy the coarse/fine parameter objects */
            MGparm_dtor(&cparm);
            MGparm_dtor(&fparm);
            Vnm_print(0, "NOsh: Done parsing ELEC section\n");
        /* Read grid dimensions */
        } else if (strcasecmp(tok, "dime") == 0) {
            /* Read the number of grid points */
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing DIME \
keyword!\n",
                  tok);
                VJMPERR2(0);
            } else {
                fparm->dime[0] = ti;
                cparm->dime[0] = ti;
            }
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing DIME \
keyword!\n",
                  tok);
                VJMPERR2(0);
            } else {
                fparm->dime[1] = ti;
                cparm->dime[1] = ti;
            }
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing DIME \
keyword!\n",
                  tok);
                VJMPERR2(0);
            } else {
                fparm->dime[2] = ti;
                cparm->dime[2] = ti;
            }
            fparm->setdime = 1;
            cparm->setdime = 1;
        /* Coarse grid length keyword */
        } else if (strcasecmp(tok, "cglen") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN keyword!\n",
                  tok);
                VJMPERR2(0);
            } else cparm->glen[0] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN keyword!\n",
                  tok);
                VJMPERR2(0);
            } else cparm->glen[1] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN keyword!\n",
                  tok);
                VJMPERR2(0);
            } else cparm->glen[2] = tf;
            cparm->setglen = 1;
        /* Fine grid length keyword */
        } else if (strcasecmp(tok, "fglen") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN
keyword!\n",
                  tok);
                VJMPERR2(0);
            } else fparm->glen[0] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN
keyword!\n",
                  tok);
                VJMPERR2(0);
            } else fparm->glen[1] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN
keyword!\n",
                  tok);
                VJMPERR2(0);
            } else fparm->glen[2] = tf;
            fparm->setglen = 1;
        /* Coarse grid center keyword */
        } else if (strcasecmp(tok, "cgcent") == 0) {
            /* If the next token isn't a float, it probably means we want to
             * center on a molecule */
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                if (strcasecmp(tok, "mol") == 0) {
                    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                    if (sscanf(tok, "%d", &ti) == 0) {
                        Vnm_print(2, "NOsh:  Read non-int (%s) while parsing GCENT MOL keyword!\n",
                          tok);
                        VJMPERR2(0);
                    } else {
                        cparm->cmeth = 1;
                        cparm->centmol = ti;
                    }
                } else {
                    Vnm_print(2, "NOsh:  Unexpected keyword (%s) while parsing GCENT!\n",
                      tok);
                    VJMPERR2(0);
                }
            } else {
                cparm->center[0] = tf;
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%lf", &tf) == 0) {
                    Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GCENT keyword!\n",
                      tok);
                    VJMPERR2(0);
                } 
                cparm->center[1] = tf;
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%lf", &tf) == 0) {
                    Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GCENT keyword!\n",
                      tok);
                    VJMPERR2(0);
                } 
                cparm->center[2] = tf;
            }   
            cparm->setgcent = 1;
        /* Fine grid center keyword */
        } else if (strcasecmp(tok, "fgcent") == 0) {
            /* If the next token isn't a float, it probably means we want to
             * center on a molecule */
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                if (strcasecmp(tok, "mol") == 0) {
                    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                    if (sscanf(tok, "%d", &ti) == 0) {
                        Vnm_print(2, "NOsh:  Read non-int (%s) while parsing GCENT MOL keyword!\n",
                          tok);
                        VJMPERR2(0);
                    } else {
                        fparm->cmeth = 1;
                        fparm->centmol = ti;
                    }
                } else {
                    Vnm_print(2, "NOsh:  Unexpected keyword (%s) while parsing GCENT!\n",
                      tok);
                    VJMPERR2(0);
                }
            } else {
                fparm->center[0] = tf;
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%lf", &tf) == 0) {
                    Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GCENT keyword!\n",
                      tok);
                    VJMPERR2(0);
                } 
                fparm->center[1] = tf;
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%lf", &tf) == 0) {
                    Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GCENT keyword!\n",
                      tok);
                    VJMPERR2(0);
                } 
                fparm->center[2] = tf;
            }   
            fparm->setgcent = 1;
        /* Read mol ID */
        } else if (strcasecmp(tok, "mol") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-int (%s) while parsing MOL keyword!\n",
                  tok);
                VJMPERR2(0);
            } 
            fparm->molid = ti;
            fparm->setmolid = 1;
            cparm->molid = ti;
            cparm->setmolid = 1;
        /* Linearized vs. nonlinear PBE */
        } else if (strcasecmp(tok, "lpbe") == 0) {
            fparm->nonlin = 0;
            fparm->setnonlin = 1;
            cparm->nonlin = 0;
            cparm->setnonlin = 1;
        } else if (strcasecmp(tok, "npbe") == 0) {
            fparm->nonlin = 1;
            fparm->setnonlin = 1;
            cparm->nonlin = 1;
            cparm->setnonlin = 1;
        /* Boundary condition flag */
        } else if (strcasecmp(tok, "bcfl") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-int (%s) while parsing BCFL keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->bcfl = ti;
            fparm->setbcfl = 1;
            cparm->bcfl = ti;
            cparm->setbcfl = 1;
        /* Ions */
        } else if (strcasecmp(tok, "ion") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->ionq[fparm->nion] = tf;
            cparm->ionq[cparm->nion] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->ionc[fparm->nion] = tf;
            cparm->ionc[cparm->nion] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->ionr[fparm->nion] = tf;
            cparm->ionr[cparm->nion] = tf;
            fparm->setion[fparm->nion] = 1;
            cparm->setion[cparm->nion] = 1;
            (fparm->nion)++;
            (cparm->nion)++;
            fparm->setnion = 1;
            cparm->setnion = 1;
        /* Solute dielectric */
        } else if (strcasecmp(tok, "pdie") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing PDIE keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->pdie = tf;
            fparm->setpdie = 1;
            cparm->pdie = tf;
            cparm->setpdie = 1;
        /* Solvent dielectric */
        } else if (strcasecmp(tok, "sdie") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SDIE keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->sdie = tf;
            fparm->setsdie = 1;
            cparm->sdie = tf;
            cparm->setsdie = 1;
        /* Surface definition methods */
        } else if (strcasecmp(tok, "srfm") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-int (%s) while parsing SRFM keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->srfm = ti;
            fparm->setsrfm = 1;
            cparm->srfm = ti;
            cparm->setsrfm = 1;
        /* Solvent radius */
        } else if (strcasecmp(tok, "srad") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SRAD keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->srad = tf;
            fparm->setsrad = 1;
            cparm->srad = tf;
            cparm->setsrad = 1;
        /* Spline window */
        } else if (strcasecmp(tok, "swin") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SWIN keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->swin = tf;
            fparm->setswin = 1;
            cparm->swin = tf;
            cparm->setswin = 1;
        /* Temperature */
        } else if (strcasecmp(tok, "temp") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing TEMP keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->temp = tf;
            fparm->settemp = 1;
            cparm->temp = tf;
            cparm->settemp = 1;
        /* Surface tension (apolar) */
        } else if (strcasecmp(tok, "gamma") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GAMMA keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->gamma = tf;
            fparm->setgamma = 1;
            cparm->gamma = tf;
            cparm->setgamma = 1;
        /* Energy writing */
        } else if (strcasecmp(tok, "calcenergy") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing WRITEENERGY keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->calcenergy = ti;
            fparm->setcalcenergy = 1;
            cparm->calcenergy = ti;
            cparm->setcalcenergy = 1;
        /* Force writing */
        } else if (strcasecmp(tok, "calcforce") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing WRITEFORCE keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->calcforce = ti;
            fparm->setcalcforce = 1;
            cparm->calcforce = ti;
            cparm->setcalcforce = 1;
        /* Potential writing */
        } else if (strcasecmp(tok, "writepot") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing WRITEPOT keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->writepot = ti;
            cparm->writepot = ti;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (strcasecmp(tok, "dx") == 0) {
                fparm->writepotfmt = 0;
                cparm->writepotfmt = 0;
            } else if (strcasecmp(tok, "avs") == 0) {
                fparm->writepotfmt = 1;
                cparm->writepotfmt = 1;
            } else if (strcasecmp(tok, "uhbd") == 0) {
                fparm->writepotfmt = 2;
                cparm->writepotfmt = 2;
            } else {
                Vnm_print(2, "NOsh:  Invalid format (%s) while parsing WRITEPOT keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            strncpy(cparm->writepotstem, tok, VMAX_ARGLEN); 
            strncpy(fparm->writepotstem, tok, VMAX_ARGLEN); 
            fparm->setwritepot = 1;
            cparm->setwritepot = 1;
         /* Accessibility writing */
        } else if (strcasecmp(tok, "writeacc") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing WRITEACC keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            fparm->writeacc = ti;
            cparm->writeacc = ti;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (strcasecmp(tok, "dx") == 0) {
                cparm->writeaccfmt = 0;
                fparm->writeaccfmt = 0;
            } else if (strcasecmp(tok, "avs") == 0) {
                fparm->writeaccfmt = 1;
                cparm->writeaccfmt = 1;
            } else if (strcasecmp(tok, "uhbd") == 0) {
                fparm->writeaccfmt = 2;
                cparm->writeaccfmt = 2;
            } else {
                Vnm_print(2, "NOsh:  Invalid format (%s) while parsing \
WRITEPOT keyword!\n",
                  tok);
                VJMPERR2(0);
            }
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            strncpy(fparm->writeaccstem, tok, VMAX_ARGLEN);
            strncpy(cparm->writeaccstem, tok, VMAX_ARGLEN);
            fparm->setwriteacc = 1;
            cparm->setwriteacc = 1;
        } else {
            Vnm_print(2, "NOsh:  Ignoring unknown keyword (%s) while parsing ELEC section!\n",
              tok);
        }
    }
  
    /* Ran out of tokens? */
    Vnm_print(2, "NOsh:  Ran out of tokens while parsing ELEC section!\n");
    VJMPERR2(0);

    VERROR1:
       Vnm_print(2, "NOsh:  Ran out of tokens while parsing ELEC section!\n");
       VJMPERR2(0);

    /* The default error handler (manages garbage collection on failure) */
    VERROR2:
        MGparm_dtor(&cparm);
        MGparm_dtor(&fparm);
        return 0;
} 
