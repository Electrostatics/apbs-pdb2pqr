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

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Private method declaration
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int NOsh_parseREAD(NOsh *thee, Vio *sock);
VPRIVATE int NOsh_parsePRINT(NOsh *thee, Vio *sock);
VPRIVATE int NOsh_parseELEC(NOsh *thee, Vio *sock);
VPRIVATE int NOsh_parseFEM(NOsh *thee, Vio *sock, NOsh_femparm *parm);
VPRIVATE int NOsh_parseMG(NOsh *thee, Vio *sock, NOsh_mgparm *parm);

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

    int i,j;

    if (thee == VNULL) return 0;

    thee->parsed = 0;
    thee->ncalc = 0;
    thee->imgcalc = 0;
    thee->ifemcalc = 0;
    thee->nmgcalc = 0;
    thee->nfemcalc = 0;
    thee->nmol = 0;
    thee->nprint = 0;

    for (j=0; j<NOSH_MAXMGPARM; j++) {
        thee->mgparm[j].setdime = 0;
        thee->mgparm[j].setnlev = 0;
        thee->mgparm[j].setgrid = 0;
        thee->mgparm[j].setglen = 0;
        thee->mgparm[j].setgcent = 0;  
        thee->mgparm[j].setmolid = 0;
        thee->mgparm[j].setnonlin = 0;
        thee->mgparm[j].setbcfl = 0;
        thee->mgparm[j].setnion = 0;
        for (i=0; i<MAXION; i++) thee->mgparm[j].setion[i] = 0;
        thee->mgparm[j].setpdie = 0;
        thee->mgparm[j].setsdie = 0;
        thee->mgparm[j].setsrfm = 0;
        thee->mgparm[j].setsrad = 0;
        thee->mgparm[j].setswin = 0; 
        thee->mgparm[j].settemp = 0;
        thee->mgparm[j].setgamma = 0;
        thee->mgparm[j].setcalcenergy = 0;      
        thee->mgparm[j].setcalcforce = 0;       
        thee->mgparm[j].setwritepot = 0; 
        thee->mgparm[j].setwriteacc = 0; 
        thee->mgparm[j].nion = 0;
        thee->mgparm[j].swin = 0;
        thee->mgparm[j].srad = 1.4;
    }

    return 1; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  NOsh_dtor
//
// Purpose:  Destroy the charge-simplex map.
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
// Purpose:  Destroy the atom object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void NOsh_dtor2(NOsh *thee) { ; }

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
            Vnm_print(0, "NOsh: Done parsing READ section\n");
        } else if (strcasecmp(tok, "print") == 0) {
            Vnm_print(0, "NOsh: Parsing PRINT section\n");
            if (!NOsh_parsePRINT(thee, sock)) return 0;
            Vnm_print(0, "NOsh: Done parsing PRINT section\n");
        } else if (strcasecmp(tok, "elec") == 0) {
            Vnm_print(0, "NOsh: Parsing ELEC section\n");
            if (!NOsh_parseELEC(thee, sock)) return 0;
            Vnm_print(0, "NOsh: Done parsing ELEC section\n");
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
                strncpy(thee->molpath[thee->nmol], tok, NOSH_MAXPATH);
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

    /* Update the calculation number */
    if (thee->ncalc >= NOSH_MAXCALC) {
        Vnm_print(2, "NOsh:  Too many electrostatics calculations in this \
run!\n");
        Vnm_print(2, "NOsh:  Current max is %d; ignoring this calculation\n", 
          NOSH_MAXCALC);
        return 1;
    }
    (thee->ncalc)++;

    if (thee->nmgcalc >= NOSH_MAXMGPARM) {
        Vnm_print(2, "NOsh:  Too many multigrid electrostatics calculations \
in this run!\n");
        Vnm_print(2, "NOsh:  Current max is %d; ignoring this calculation\n", 
          NOSH_MAXMGPARM);
        return 1;
    }

    /* The next token HAS to be the method */
    if (Vio_scanf(sock, "%s", tok) == 1) {
        /* THERE ARE TWO CHOICES HERE: MG-MANUAL OR MG-AUTOMATIC
         * IN THE CASE OF THE LATTER, SEVERAL FOCUSSING CALCULATIONS ARE
         * PERFORMED TO ACHEIVE THE TARGET RESOLUTION.  THIS MEANS THAT THINGS
         * LIKE (thee->nmgcalc) AND ARRAYS WHICH ARE REFERENCED BY THIS
         * VARIABLE CAN NO LONGER BE UPDATED OUTSIDE OF THE PARSING ROUTINE.
         * IT ALSO MEANS THAT PARSING OF THE PRINT STATEMENTS WILL HAVE TO BE
         * HANDLED CAREFULLY TO USE ONLY THE SOLUTION FROM THE FINEST LEVEL OF
         * THE AUTOMATIC FOCUSSING HIERARCHY */
        if (strcasecmp(tok, "mg") == 0) {
	    /* Check to see if he have any room left for this type of
             * calculation, if so: set the calculation type, update the number
             * of calculations of this type, and parse the rest of the section
             */
            if (thee->nmgcalc >= NOSH_MAXMGPARM) {
                Vnm_print(2, "NOsh:  Too many multigrid electrostatics \
calculations in this run!\n");
                Vnm_print(2, "NOsh:  Current max is %d; ignoring this \
calculation\n",
                  NOSH_MAXMGPARM);
                return 1;
            }
            thee->calctype[thee->ncalc - 1] = 0;
            (thee->nmgcalc)++;
            Vnm_print(0, "NOsh: Parsing parameters for MG calculation #%d\n",
              thee->nmgcalc);
            return NOsh_parseMG(thee, sock, &(thee->mgparm[thee->nmgcalc-1]));
        } else if (strcasecmp(tok, "fem") == 0) {
            /* Check to see if he have any room left for this type of
             * calculation, if so: set the calculation type, update the number
             * of calculations of this type, and parse the rest of the section
             */
            if (thee->nmgcalc >= NOSH_MAXFEMPARM) {
                Vnm_print(2, "NOsh:  Too many finite element electrostatics calculations in this run!\n");
                Vnm_print(2, "NOsh:  Current max is %d; ignoring this calculation\n",
                  NOSH_MAXFEMPARM);
                return 1;
            }
            thee->calctype[thee->ncalc - 1] = 1;
            (thee->nfemcalc)++;
            Vnm_print(0, "NOsh: Parsing parameters for FEM calculation #%d\n",
              thee->nfemcalc);
            return NOsh_parseFEM(thee,sock,&(thee->femparm[thee->nfemcalc-1]));
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
VPRIVATE int NOsh_parseFEM(NOsh *thee, Vio *sock, NOsh_femparm *parm) {

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
// Routine:  NOsh_parseMG
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
VPRIVATE int NOsh_parseMG(NOsh *thee, Vio *sock, NOsh_mgparm *parm) {

    char tok[VMAX_BUFSIZE];
    double tf;
    int i, ti;

    if (thee == VNULL) {
        Vnm_print(2, "NOsh_parseMG:  Got NULL thee!\n");
        return 0;
    }

    if (parm == VNULL) {
        Vnm_print(2, "NOsh_parseMG:  Got NULL parm!\n");
        return 0;
    }

    if (sock == VNULL) {
        Vnm_print(2, "NOsh_parseMG:  Got pointer to NULL socket!\n");
        return 0;
    }

    if (thee->parsed) {
        Vnm_print(2, "NOsh_parseMG:  Already parsed an input file!\n");
        return 0;
    }

    /* Here we go... */
    while (Vio_scanf(sock, "%s", tok) == 1) {
        if (strcasecmp(tok, "end") == 0) {
            /* Check to see that everything was set */
            if (!parm->setdime) {
                Vnm_print(2, "NOsh: DIME not set!\n");
                return 0;
            }
            if (!parm->setnlev) {
                Vnm_print(2, "NOsh: NLEV not set!\n");
                return 0;
            }
            if ((!parm->setgrid) && (!parm->setglen)) {
                Vnm_print(2, "NOsh: Neither GRID nor GLEN set!\n");
                return 0;
            }
            if ((parm->setgrid) && (parm->setglen)) {
                Vnm_print(2, "NOsh: Both GRID and GLEN set!\n");
                return 0;
            }
            if (!parm->setgcent) {
                Vnm_print(2, "NOsh: GCENT not set!\n");
                return 0;
            }
            if (!parm->setmolid) {
                Vnm_print(2, "NOsh: MOL not set!\n");
                return 0;
            }
            if (!parm->setnonlin) {
                Vnm_print(2, "NOsh: LPBE or NPBE not set!\n");
                return 0;
            }
            if (!parm->setbcfl) {
                Vnm_print(2, "NOsh: BCFL not set!\n");
                return 0;
            }
            if (!parm->setnion) {
                parm->setnion = 1;
                parm->nion = 0;
            } else {
                if (parm->nion > 2) {
                    Vnm_print(2, "NOsh:  Only 1:1 ionic solutions presently allowed!\n");           
                    return 0;
                }
                Vnm_print(2, "NOsh:  Warning -- only the largest ionic radius is used!\n");
            }
            for (i=0; i<parm->nion; i++) {
                if (!parm->setion[i]) {
                    Vnm_print(2, "NOsh: ION #%d not set!\n",i);
                    return 0;
                }
                if (VABS(parm->ionq[i]) != 1.00) {
                    Vnm_print(2, "NOsh: Only monovalent ions presently allowed!\n");
                    return 0;
                }
            }
            if (!parm->setpdie) {
                Vnm_print(2, "NOsh: PDIE not set!\n");
                return 0;
            }
            if (!parm->setsdie) {
                Vnm_print(2, "NOsh: SDIE not set!\n");
                return 0;
            }
            if (!parm->setsrfm) {
                Vnm_print(2, "NOsh: SRFM not set!\n");
                return 0;
            }
            if (((parm->srfm==0) || (parm->srfm==1)) && (!parm->setsrad)) {
                Vnm_print(2, "NOsh: SRAD not set!\n");
                return 0;
            }
            if ((parm->srfm==2) && (!parm->setswin)) {
                Vnm_print(2, "NOsh: SWIN not set!\n");
                return 0;
            }
            if (!parm->settemp) {
                Vnm_print(2, "NOsh: TEMP not set!\n");
                return 0;
            }
            if (!parm->setgamma) {
                Vnm_print(2, "NOsh: GAMMA not set!\n");
                return 0;
            }
            if (!parm->setcalcenergy) parm->calcenergy = 0;
            if (!parm->setcalcforce) parm->calcforce = 0;
            if (!parm->setwritepot) parm->writepot = 0;
            if (!parm->setwriteacc) parm->writeacc = 0;
            Vnm_print(0, "NOsh: Done parsing ELEC section\n");
            return 1;
        /* Read grid dimensions */
        } else if (strcasecmp(tok, "dime") == 0) {
            /* Read the number of grid points */
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing DIME keyword!\n",
                  tok);
                return 0;
            } else parm->dime[0] = ti;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing DIME keyword!\n",
                  tok);
                return 0;
            } else parm->dime[1] = ti;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing DIME keyword!\n",
                  tok);
                return 0;
            } else parm->dime[2] = ti;
            parm->setdime = 1;
        /* Read number of levels in hierarchy */
        } else if (strcasecmp(tok, "nlev") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-integer (%s) while parsing NLEV keyword!\n",
                  tok);
                return 0;
            } else parm->nlev = ti;
            parm->setnlev = 1;
        } else if (strcasecmp(tok, "grid") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GRID keyword!\n",
                  tok);
                return 0;
            } else parm->grid[0] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GRID keyword!\n",
                  tok);
                return 0;
            } else parm->grid[1] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GRID keyword!\n",
                  tok);
                return 0;
            } else parm->grid[2] = tf;
            parm->setgrid = 1;
        /* Grid length keyword */
        } else if (strcasecmp(tok, "glen") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN keyword!\n",
                  tok);
                return 0;
            } else parm->glen[0] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN keyword!\n",
                  tok);
                return 0;
            } else parm->glen[1] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GLEN keyword!\n",
                  tok);
                return 0;
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
                        return 0;
                    } else {
                        parm->cmeth = 1;
                        parm->centmol = ti;
                    }
                } else {
                    Vnm_print(2, "NOsh:  Unexpected keyword (%s) while parsing GCENT!\n",
                      tok);
                    return 0;
                }
            } else {
                parm->center[0] = tf;
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%lf", &tf) == 0) {
                    Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GCENT keyword!\n",
                      tok);
                    return 0;
                } 
                parm->center[1] = tf;
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%lf", &tf) == 0) {
                    Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GCENT keyword!\n",
                      tok);
                    return 0;
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
                return 0;
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
                return 0;
            }
            parm->bcfl = ti;
            parm->setbcfl = 1;
        /* Ions */
        } else if (strcasecmp(tok, "ion") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION keyword!\n",
                  tok);
                return 0;
            }
            parm->ionq[parm->nion] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION keyword!\n",
                  tok);
                return 0;
            }
            parm->ionc[parm->nion] = tf;
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION keyword!\n",
                  tok);
                return 0;
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
                return 0;
            }
            parm->pdie = tf;
            parm->setpdie = 1;
        /* Solvent dielectric */
        } else if (strcasecmp(tok, "sdie") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SDIE keyword!\n",
                  tok);
                return 0;
            }
            parm->sdie = tf;
            parm->setsdie = 1;
        /* Surface definition methods */
        } else if (strcasecmp(tok, "srfm") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-int (%s) while parsing SRFM keyword!\n",
                  tok);
                return 0;
            }
            parm->srfm = ti;
            parm->setsrfm = 1;
        /* Solvent radius */
        } else if (strcasecmp(tok, "srad") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SRAD keyword!\n",
                  tok);
                return 0;
            }
            parm->srad = tf;
            parm->setsrad = 1;
        /* Spline window */
        } else if (strcasecmp(tok, "swin") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SWIN keyword!\n",
                  tok);
                return 0;
            }
            parm->swin = tf;
            parm->setswin = 1;
        /* Temperature */
        } else if (strcasecmp(tok, "temp") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing TEMP keyword!\n",
                  tok);
                return 0;
            }
            parm->temp = tf;
            parm->settemp = 1;
        /* Surface tension (apolar) */
        } else if (strcasecmp(tok, "gamma") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &tf) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GAMMA keyword!\n",
                  tok);
                return 0;
            }
            parm->gamma = tf;
            parm->setgamma = 1;
        /* Energy writing */
        } else if (strcasecmp(tok, "calcenergy") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing WRITEENERGY keyword!\n",
                  tok);
                return 0;
            }
            parm->calcenergy = ti;
            parm->setcalcenergy = 1;
        /* Force writing */
        } else if (strcasecmp(tok, "calcforce") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing WRITEFORCE keyword!\n",
                  tok);
                return 0;
            }
            parm->calcforce = ti;
            parm->setcalcforce = 1;
        /* Potential writing */
        } else if (strcasecmp(tok, "writepot") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing WRITEPOT keyword!\n",
                  tok);
                return 0;
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
                return 0;
            }
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            strncpy(parm->writepotstem, tok, NOSH_MAXPATH); 
            parm->setwritepot = 1;
         /* Accessibility writing */
        } else if (strcasecmp(tok, "writeacc") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &ti) == 0) {
                Vnm_print(2, "NOsh:  Read non-float (%s) while parsing WRITEACC keyword!\n",
                  tok);
                return 0;
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
                return 0;
            }
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            strncpy(parm->writeaccstem, tok, NOSH_MAXPATH);
            parm->setwriteacc = 1;
        } else {
            Vnm_print(2, "NOsh:  Ignoring unknown keyword (%s) while parsing ELEC section!\n",
              tok);
        }
    }
  
    /* Ran out of tokens? */
    Vnm_print(2, "NOsh_parseMG:  Ran out of tokens while parsing ELEC section!\n");
    return 0;

    VERROR1:
       Vnm_print(2, "NOsh_parseMG:  Ran out of tokens while parsing ELEC \
section!\n");
       return 0;
}
