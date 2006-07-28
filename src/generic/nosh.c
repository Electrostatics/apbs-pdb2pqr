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
 * Center for Computational Biology
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2002-2006.  Washington University in St. Louis.
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
#include "apbs/nosh.h"
#include "apbs/vstring.h"
#include "apbs/mgparm.h"
#include "apbs/femparm.h"

VEMBED(rcsid="$Id$")


VPRIVATE int NOsh_parseREAD(
							NOsh *thee, 
							Vio *sock);

VPRIVATE int NOsh_parsePRINT(
							 NOsh *thee, 
							 Vio *sock);

VPRIVATE int NOsh_parseELEC(
							NOsh *thee, 
							Vio *sock
							);

VEXTERNC int NOsh_parseFEM(
						   NOsh *thee, 
						   Vio *sock, 
						   NOsh_calc *elec
						   );

VEXTERNC int NOsh_parseMG(
						  NOsh *thee, 
						  Vio *sock, 
						  NOsh_calc *elec
						  );

VPRIVATE int NOsh_setupCalcMG(
							  NOsh *thee,
							  NOsh_calc *elec
							  );	

VPRIVATE int NOsh_setupCalcMGAUTO(
								  NOsh *thee,
								  NOsh_calc *elec
								  );	

VPRIVATE int NOsh_setupCalcMGMANUAL(
									NOsh *thee,
									NOsh_calc *elec
									);	

VPRIVATE int NOsh_setupCalcMGPARA(
								  NOsh *thee,
								  NOsh_calc *elec
								  );	

VPRIVATE int NOsh_setupCalcFEM(
							   NOsh *thee,
							   NOsh_calc *elec
							   );

VPRIVATE int NOsh_setupCalcFEMANUAL(
							   NOsh *thee,
							   NOsh_calc *elec
							   );

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
	return thee->calc[icalc];
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

VPUBLIC NOsh_PrintType NOsh_printWhat(NOsh *thee, int iprint) {
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

VPUBLIC char* NOsh_elecname(NOsh *thee, int ielec) {
	VASSERT(thee != VNULL);
	VASSERT(ielec < thee->nelec + 1);
	return thee->elecname[ielec];
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

VPUBLIC NOsh* NOsh_ctor(int rank, int size) {
	
	/* Set up the structure */
	NOsh *thee = VNULL;
	thee = Vmem_malloc(VNULL, 1, sizeof(NOsh) );
	VASSERT( thee != VNULL);
	VASSERT( NOsh_ctor2(thee, rank, size) );
	
	return thee;
}

VPUBLIC int NOsh_ctor2(NOsh *thee, int rank, int size) {
	
	int i;
	
	if (thee == VNULL) return 0;
	
	thee->proc_rank = rank;
	thee->proc_size = size;
	
	thee->ispara = 0;
    thee->parsed = 0;
	
    thee->nmol = 0;
    thee->gotparm = 0;
    thee->ncharge = 0;
    thee->ndiel = 0;
    thee->nkappa = 0;
    thee->nprint = 0;
	
    for (i=0; i<NOSH_MAXCALC; i++) {
		thee->calc[i] = VNULL;
		thee->elec[i] = VNULL;
	}
	for (i=0; i<NOSH_MAXMOL; i++) {
		thee->alist[i] = VNULL;
	}
	thee->ncalc = 0;
	thee->nelec = 0;
	
    return 1; 
}

VPUBLIC void NOsh_dtor(NOsh **thee) {
    if ((*thee) != VNULL) {
        NOsh_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(NOsh), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void NOsh_dtor2(NOsh *thee) { 
	
	int i;
	
	if (thee != VNULL) {
		for (i=0; i<(thee->ncalc); i++) NOsh_calc_dtor(&(thee->calc[i]));
		for (i=0; i<(thee->nelec); i++) NOsh_calc_dtor(&(thee->elec[i]));
    }
	
}

VPUBLIC NOsh_calc* NOsh_calc_ctor(
								  NOsh_CalcType calctype
								  ) {
	NOsh_calc *thee;
	thee = (NOsh_calc *)Vmem_malloc(VNULL, 1, sizeof(NOsh_calc));
	thee->calctype = calctype;
	switch (calctype) {
		case NCT_MG:
			thee->mgparm = MGparm_ctor(MCT_NONE);
			thee->femparm = VNULL;
			break;
		case NCT_FEM:
			thee->mgparm = VNULL;
			thee->femparm = FEMparm_ctor(FCT_NONE);
			break;
		default:
			Vnm_print(2, "NOsh_calc_ctor:  unknown calculation type (%d)!\n",
					  calctype);
			VASSERT(0);
	}
	thee->pbeparm = PBEparm_ctor();
	
	return thee;
}

VPUBLIC void NOsh_calc_dtor(
							NOsh_calc **thee 
							) {
	
	NOsh_calc *calc = VNULL;
	calc = *thee;
	if (calc == VNULL) return;
	
	switch (calc->calctype) {
		case NCT_MG:
			MGparm_dtor(&(calc->mgparm));
			break;
		case NCT_FEM:
			FEMparm_dtor(&(calc->femparm));
			break;
		default:
			Vnm_print(2, "NOsh_calc_ctor:  unknown calculation type (%d)!\n",
					  calc->calctype);
			VASSERT(0);
	}
	PBEparm_dtor(&(calc->pbeparm));
	
	Vmem_free(VNULL, 1, sizeof(NOsh_calc), (void **)thee);
	calc = VNULL;
	
}

VPUBLIC int NOsh_calc_copy(
						   NOsh_calc *thee,
						   NOsh_calc *source
						   ) {
	
	VASSERT(thee != VNULL);
	VASSERT(source != VNULL);
	VASSERT(thee->calctype == source->calctype);
	if (source->mgparm != VNULL)
		MGparm_copy(thee->mgparm, source->mgparm);
	if (source->femparm != VNULL)
		FEMparm_copy(thee->femparm, source->femparm);
	if (source->pbeparm != VNULL)
		PBEparm_copy(thee->pbeparm, source->pbeparm);
	
	return 1;
	
}

VPUBLIC int NOsh_parseInputFile(
								NOsh *thee, 
								char *filename
								) {
	
	Vio *sock;
	int rc;
	
	sock = Vio_ctor("FILE", "ASC", VNULL, filename, "r");
	rc = NOsh_parseInput(thee, sock);
	Vio_dtor(&sock);
	
	return rc;
}

VPUBLIC int NOsh_parseInput(
							NOsh *thee, 
							Vio *sock
							) {
	
	char *MCwhiteChars = " =,;\t\n";
	char *MCcommChars  = "#%";
	char tok[VMAX_BUFSIZE];
	
	if (thee == VNULL) {
		Vnm_print(2, "NOsh_parseInput:  Got NULL thee!\n");
		return 0;
	}
	
	if (sock == VNULL) {
		Vnm_print(2, "NOsh_parseInput:  Got pointer to NULL socket!\n");
		Vnm_print(2, "NOsh_parseInput:  The specified input file was not found!\n");
		return 0;
	} 
	
	if (thee->parsed) {
		Vnm_print(2, "NOsh_parseInput:  Already parsed an input file!\n");
		return 0;
	}
	
	if (Vio_accept(sock, 0) < 0) {
		Vnm_print(2, "NOsh_parseInput:  Problem reading from socket!\n");
		return 0;
	}
	
	/* Set up the whitespace and comment character definitions */
	Vio_setWhiteChars(sock, MCwhiteChars);
	Vio_setCommChars(sock, MCcommChars);
	
	/* We parse the file until we run out of tokens */
	Vnm_print(0, "NOsh_parseInput:  Starting file parsing...\n");
	while (Vio_scanf(sock, "%s", tok) == 1) {
		/* At the highest level, we look for keywords that indicate functions 
like:
		
		read => Read in a molecule file
		elec => Do an electrostatics calculation
		print => Print some results
		quit => Quit
		
		These cause the code to go to a lower-level parser routine which
		handles keywords specific to the particular function.  Each
		lower-level parser routine then returns when it hits the "end"
		keyword.  Due to this simple layout, no nesting of these "function"
		sections is allowed.
		*/
		if (Vstring_strcasecmp(tok, "read") == 0) {
			Vnm_print(0, "NOsh: Parsing READ section\n");
			if (!NOsh_parseREAD(thee, sock)) return 0;
			Vnm_print(0, "NOsh: Done parsing READ section \
(nmol=%d, ndiel=%d, nkappa=%d, ncharge=%d)\n", thee->nmol, thee->ndiel, 
					  thee->nkappa, thee->ncharge);
		} else if (Vstring_strcasecmp(tok, "print") == 0) {
			Vnm_print(0, "NOsh: Parsing PRINT section\n");
			if (!NOsh_parsePRINT(thee, sock)) return 0;
			Vnm_print(0, "NOsh: Done parsing PRINT section\n");
		} else if (Vstring_strcasecmp(tok, "elec") == 0) {
			Vnm_print(0, "NOsh: Parsing ELEC section\n");
			if (!NOsh_parseELEC(thee, sock)) return 0;
			Vnm_print(0, "NOsh: Done parsing ELEC section (nelec = %d)\n",
					  thee->nelec);
		} else if (Vstring_strcasecmp(tok, "quit") == 0) {
			Vnm_print(0, "NOsh: Done parsing file (got QUIT)\n");
			break;
		} else {
			Vnm_print(2, "NOsh_parseInput: Ignoring undefined keyword %s!\n", tok);
		}
	}
	
	thee->parsed = 1;
	return 1;
	
}

VPRIVATE int NOsh_parseREAD_MOL(NOsh *thee, Vio *sock) {
	
    char tok[VMAX_BUFSIZE];
    NOsh_MolFormat molfmt;
	
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (Vstring_strcasecmp(tok, "pqr") == 0) {
        molfmt = NMF_PQR;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        Vnm_print(0, "NOsh: Storing molecule %d path %s\n", 
				  thee->nmol, tok);
        thee->molfmt[thee->nmol] = molfmt;
        strncpy(thee->molpath[thee->nmol], tok, VMAX_ARGLEN);
        (thee->nmol)++;
    } else if (Vstring_strcasecmp(tok, "pdb") == 0) {
        molfmt = NMF_PDB;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        Vnm_print(0, "NOsh: Storing molecule %d path %s\n", 
				  thee->nmol, tok);
        thee->molfmt[thee->nmol] = molfmt;
        strncpy(thee->molpath[thee->nmol], tok, VMAX_ARGLEN);
        (thee->nmol)++;
    } else if (Vstring_strcasecmp(tok, "xml") == 0) {
        molfmt = NMF_XML;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        Vnm_print(0, "NOsh: Storing molecule %d path %s\n", 
				  thee->nmol, tok);
        thee->molfmt[thee->nmol] = molfmt;
        strncpy(thee->molpath[thee->nmol], tok, VMAX_ARGLEN);
        (thee->nmol)++;
    } else {
        Vnm_print(2, "NOsh_parseREAD:  Ignoring undefined mol format \
%s!\n", tok);
    } 
	
    return 1;
	
	
VERROR1:
        Vnm_print(2, "NOsh_parseREAD_MOL:  Ran out of tokens while parsing READ section!\n");
	return 0;
	
}

VPRIVATE int NOsh_parseREAD_PARM(NOsh *thee, Vio *sock) {
	
    char tok[VMAX_BUFSIZE];
    NOsh_ParmFormat parmfmt;
	
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (Vstring_strcasecmp(tok, "flat") == 0) {
        parmfmt = NPF_FLAT;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (thee->gotparm) {
            Vnm_print(2, "NOsh:  Hey!  You already specified a parameterfile (%s)!\n", thee->parmpath);
            Vnm_print(2, "NOsh:  I'm going to ignore this one (%s)!\n", tok);
        } else {
            thee->parmfmt = parmfmt;
            thee->gotparm = 1;
            strncpy(thee->parmpath, tok, VMAX_ARGLEN);
        }
    } else if(Vstring_strcasecmp(tok, "xml") == 0) {
        parmfmt = NPF_XML;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (thee->gotparm) {
            Vnm_print(2, "NOsh:  Hey!  You already specified a parameterfile (%s)!\n", thee->parmpath);
            Vnm_print(2, "NOsh:  I'm going to ignore this one (%s)!\n", tok);
        } else {
            thee->parmfmt = parmfmt;
            thee->gotparm = 1;
            strncpy(thee->parmpath, tok, VMAX_ARGLEN);
        }
		
    } else {
        Vnm_print(2, "NOsh_parseREAD:  Ignoring undefined parm format \
%s!\n", tok);
    } 
	
    return 1;
	
VERROR1:
        Vnm_print(2, "NOsh_parseREAD_PARM:  Ran out of tokens while parsing READ section!\n");
	return 0;
	
}

VPRIVATE int NOsh_parseREAD_DIEL(NOsh *thee, Vio *sock) {
	
    char tok[VMAX_BUFSIZE];
    Vdata_Format dielfmt;
	
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (Vstring_strcasecmp(tok, "dx") == 0) {
        dielfmt = VDF_DX;
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
	
    return 1;
	
VERROR1:
        Vnm_print(2, "NOsh_parseREAD_DIEL:  Ran out of tokens while parsing READ \
section!\n");
	return 0;
	
}

VPRIVATE int NOsh_parseREAD_KAPPA(NOsh *thee, Vio *sock) {
	
    char tok[VMAX_BUFSIZE];
    Vdata_Format kappafmt;
	
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (Vstring_strcasecmp(tok, "dx") == 0) {
        kappafmt = VDF_DX;
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
	
    return 1;
	
VERROR1:
        Vnm_print(2, "NOsh_parseREAD:  Ran out of tokens while parsing READ \
section!\n");
	return 0;
	
}

VPRIVATE int NOsh_parseREAD_CHARGE(NOsh *thee, Vio *sock) {
	
    char tok[VMAX_BUFSIZE];
    Vdata_Format chargefmt;
	
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (Vstring_strcasecmp(tok, "dx") == 0) {
        chargefmt = VDF_DX;
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
	
    return 1;
	
VERROR1:
        Vnm_print(2, "NOsh_parseREAD:  Ran out of tokens while parsing READ \
section!\n");
	return 0;
	
}

VPRIVATE int NOsh_parseREAD(NOsh *thee, Vio *sock) {
	
    char tok[VMAX_BUFSIZE];
	
    if (thee == VNULL) {
        Vnm_print(2, "NOsh_parseREAD:  Got NULL thee!\n");
        return 0;
    }
	
    if (sock == VNULL) {
        Vnm_print(2, "NOsh_parseREAD:  Got pointer to NULL socket!\n");
        return 0;
    } 
	
    if (thee->parsed) {
        Vnm_print(2, "NOsh_parseREAD:  Already parsed an input file!\n");
        return 0;
    }
	
    /* Read until we run out of tokens (bad) or hit the "END" keyword (good) */
    while (Vio_scanf(sock, "%s", tok) == 1) {
		
        if (Vstring_strcasecmp(tok, "end") == 0) {
            Vnm_print(0, "NOsh: Done parsing READ section\n");
            return 1;
        } else if (Vstring_strcasecmp(tok, "mol") == 0) {
            NOsh_parseREAD_MOL(thee, sock);
        } else if (Vstring_strcasecmp(tok, "parm") == 0) {
            NOsh_parseREAD_PARM(thee,sock);
        } else if (Vstring_strcasecmp(tok, "diel") == 0) {
            NOsh_parseREAD_DIEL(thee,sock);
        } else if (Vstring_strcasecmp(tok, "kappa") == 0) {
            NOsh_parseREAD_KAPPA(thee,sock);
        } else if (Vstring_strcasecmp(tok, "charge") == 0) {
            NOsh_parseREAD_CHARGE(thee,sock);
        } else {
            Vnm_print(2, "NOsh_parseREAD:  Ignoring undefined keyword %s!\n", 
					  tok);
        }
    }
	
	
    /* We ran out of tokens! */
    Vnm_print(2, "NOsh_parseREAD:  Ran out of tokens while parsing READ \
section!\n");
    return 0;
	
}

VPRIVATE int NOsh_parsePRINT(NOsh *thee, Vio *sock) {
	
    char tok[VMAX_BUFSIZE];
	char name[VMAX_BUFSIZE];
    int ti, idx, expect, ielec;
	
    if (thee == VNULL) {
        Vnm_print(2, "NOsh_parsePRINT:  Got NULL thee!\n");
        return 0;
    }
	
    if (sock == VNULL) {
        Vnm_print(2, "NOsh_parsePRINT:  Got pointer to NULL socket!\n");
        return 0;
    } 
	
    if (thee->parsed) {
        Vnm_print(2, "NOsh_parsePRINT:  Already parsed an input file!\n");
        return 0;
    }
	
    idx = thee->nprint;
    if (thee->nprint >= NOSH_MAXPRINT) {
        Vnm_print(2, "NOsh_parsePRINT:  Exceeded max number (%d) of PRINT \
sections\n", 
				  NOSH_MAXPRINT);
        return 0;
    }
	
	
    /* The first thing we read is the thing we want to print */ 
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (Vstring_strcasecmp(tok, "energy") == 0) {
        thee->printwhat[idx] = NPT_ENERGY;
        thee->printnarg[idx] = 0;
    } else if (Vstring_strcasecmp(tok, "force") == 0) {
        thee->printwhat[idx] = NPT_FORCE;
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
            if ((sscanf(tok, "%d", &ti) == 1) && 
				(Vstring_isdigit(tok) == 1)) {
                if (expect == 0) {
                    thee->printcalc[idx][thee->printnarg[idx]] = ti-1;
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
				/* Grab a calculation name from elec ID */
			} else if (sscanf(tok, "%s", name) == 1) {
			    if (expect == 0) {
					for (ielec=0; ielec<thee->nelec; ielec++) {
					    if (Vstring_strcasecmp(thee->elecname[ielec], name) == 0) {
							thee->printcalc[idx][thee->printnarg[idx]] = ielec;
							expect = 1;
							break;
					    }
				    }
					if (expect == 0) {
						Vnm_print(2, "No ELEC statement has been named %s!\n", 
								  name);
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

VPRIVATE int NOsh_parseELEC(NOsh *thee, Vio *sock) {
	
	NOsh_calc *calc = VNULL;
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
	
    /* Get a pointer to the latest ELEC calc object and update the ELEC 
		statement number */
    if (thee->nelec >= NOSH_MAXCALC) {
        Vnm_print(2, "NOsh:  Too many electrostatics calculations in this \
run!\n");
        Vnm_print(2, "NOsh:  Current max is %d; ignoring this calculation\n", 
				  NOSH_MAXCALC);
        return 1;
    }
	
	
	
    /* The next token HAS to be the method OR "name" */
	if (Vio_scanf(sock, "%s", tok) == 1) {
	    if (Vstring_strcasecmp(tok, "name") == 0) {
		    Vio_scanf(sock, "%s", tok);
            strncpy(thee->elecname[thee->nelec], tok, VMAX_ARGLEN);
			if (Vio_scanf(sock, "%s", tok) != 1) {
				Vnm_print(2, "NOsh_parseELEC:  Ran out of tokens while reading \
ELEC section!\n");
				return 0;
			}
		}
		if (Vstring_strcasecmp(tok, "mg-manual") == 0) {
			thee->elec[thee->nelec] = NOsh_calc_ctor(NCT_MG);
			calc = thee->elec[thee->nelec];
			(thee->nelec)++;	
			calc->mgparm->type = MCT_MANUAL;
		    return NOsh_parseMG(thee, sock, calc);
		} else if (Vstring_strcasecmp(tok, "mg-auto") == 0) {
			thee->elec[thee->nelec] = NOsh_calc_ctor(NCT_MG);
			calc = thee->elec[thee->nelec];
			(thee->nelec)++;	
			calc->mgparm->type = MCT_AUTO;
		    return NOsh_parseMG(thee, sock, calc);
		} else if (Vstring_strcasecmp(tok, "mg-para") == 0) {
			thee->elec[thee->nelec] = NOsh_calc_ctor(NCT_MG);
			calc = thee->elec[thee->nelec];
			(thee->nelec)++;	
			calc->mgparm->type = MCT_PARALLEL;
		    return NOsh_parseMG(thee, sock, calc);
		} else if (Vstring_strcasecmp(tok, "mg-dummy") == 0) {
			thee->elec[thee->nelec] = NOsh_calc_ctor(NCT_MG);
			calc = thee->elec[thee->nelec];
			(thee->nelec)++;	
			calc->mgparm->type = MCT_DUMMY;
		    return NOsh_parseMG(thee, sock, calc);
		} else if (Vstring_strcasecmp(tok, "fe-manual") == 0) {
			thee->elec[thee->nelec] = NOsh_calc_ctor(NCT_FEM);
			calc = thee->elec[thee->nelec];
			(thee->nelec)++;	
			calc->femparm->type = FCT_MANUAL;
		    return NOsh_parseFEM(thee, sock, calc);
		} else {
		    Vnm_print(2, "NOsh_parseELEC: The method (\"mg\" or \"fem\") or \
\"name\" must be the first keyword in the ELEC section\n");
			return 0;
		} 
	} 
	
    Vnm_print(2, "NOsh_parseELEC:  Ran out of tokens while reading ELEC section!\n");
    return 0;
	
}

VPUBLIC int NOsh_setupCalc(
						   NOsh *thee, 
						   Valist *alist[NOSH_MAXMOL]
						   ) {
	int ielec, imol;
	NOsh_calc *elec = VNULL;
	
	VASSERT(thee != VNULL);
	for (imol=0; imol<thee->nmol; imol++) {
		thee->alist[imol] = alist[imol];
	}
	
	
	for (ielec=0; ielec<(thee->nelec); ielec++) {
		/* Unload the calculation object containing the ELEC information */
		elec = thee->elec[ielec];
		
		
		/* Check the calculation type */
		switch (elec->calctype) {
			case NCT_MG:
				NOsh_setupCalcMG(thee, elec);
				break;
			case NCT_FEM:
				NOsh_setupCalcFEM(thee, elec);
				break;
			default:
				Vnm_print(2, "NOsh_setupCalc:  Invalid calculation type (%d)!\n",
						  elec->calctype);
				return 0;
		}
		
		/* At this point, the most recently-created NOsh_calc object should be the
			one we use for results for this ELEC statement.  Assign it. */
		/* Associate ELEC statement with the calculation */
		thee->elec2calc[ielec] = thee->ncalc-1;
		Vnm_print(0, "NOsh_setupCalc:  Mapping ELEC statement %d (%d) to \
calculation %d (%d)\n", ielec, ielec+1, thee->elec2calc[ielec], 
				  thee->elec2calc[ielec]+1);
	}
	
	return 1;
}



VPUBLIC int NOsh_parseMG(
						 NOsh *thee, 
						 Vio *sock, 
						 NOsh_calc *elec
						 ) {
	
	char tok[VMAX_BUFSIZE];
	MGparm *mgparm = VNULL;    
	PBEparm *pbeparm = VNULL;    
	int rc;
    
	/* Check the arguments */
	if (thee == VNULL) {
		Vnm_print(2, "NOsh:  Got NULL thee!\n");
		return 0;
	}
	if (sock == VNULL) {
		Vnm_print(2, "NOsh:  Got pointer to NULL socket!\n");
		return 0;
	}
	if (elec == VNULL) {
		Vnm_print(2, "NOsh:  Got pointer to NULL elec object!\n");
		return 0;
	}
	mgparm = elec->mgparm;
	if (mgparm == VNULL) {
		Vnm_print(2, "NOsh:  Got pointer to NULL mgparm object!\n");
		return 0;
	}
	pbeparm = elec->pbeparm;
	if (pbeparm == VNULL) {
		Vnm_print(2, "NOsh:  Got pointer to NULL pbeparm object!\n");
		return 0;
	}
	
	Vnm_print(0, "NOsh_parseMG: Parsing parameters for MG calculation\n");
	
	/* Parallel stuff */
	if (mgparm->type == MCT_PARALLEL) {
		mgparm->proc_rank = thee->proc_rank;
		mgparm->proc_size = thee->proc_size;
		mgparm->setrank = 1;
		mgparm->setsize = 1;
	}

	
	/* Start snarfing tokens from the input stream */
	rc = 1;
	while (Vio_scanf(sock, "%s", tok) == 1) { 
		
		Vnm_print(0, "NOsh_parseMG:  Parsing %s...\n", tok);
		
		/* See if it's an END token */
		if (Vstring_strcasecmp(tok, "end") == 0) {
			mgparm->parsed = 1;
			pbeparm->parsed = 1;
			rc = 1;
			break;
		}
		
		/* Pass the token through a series of parsers */
		rc = PBEparm_parseToken(pbeparm, tok, sock);
		if (rc == -1) {
			Vnm_print(0, "NOsh_parseMG:  parsePBE error!\n");
			break;
		} else if (rc == 0) {
			/* Pass the token to the generic MG parser */
			rc = MGparm_parseToken(mgparm, tok, sock);
			if (rc == -1) { 
				Vnm_print(0, "NOsh_parseMG:  parseMG error!\n");
				break;
			} else if (rc == 0) {
				/* We ran out of parsers! */
				Vnm_print(2, "NOsh:  Unrecognized keyword: %s\n", tok);
				break;
			}
		}
	}
	
	/* Handle various errors arising in the token-snarfing loop -- these all 
		just result in simple returns right now */
	if (rc == -1) return 0;
	if (rc == 0) return 0;
	
	/* Check the status of the parameter objects */
	if ((!MGparm_check(mgparm)) || (!PBEparm_check(pbeparm))) {
		Vnm_print(2, "NOsh:  MG parameters not set correctly!\n");
		return 0;
	}
	
	return 1;
}

VPRIVATE int NOsh_setupCalcMG(
							  NOsh *thee,
							  NOsh_calc *calc
							  ) {
	
	MGparm *mgparm = VNULL;
	
	VASSERT(thee != VNULL);
	VASSERT(calc != VNULL);
	mgparm = calc->mgparm;
	VASSERT(mgparm != VNULL);
	
	
	/* Now we're ready to whatever sorts of post-processing operations that are
		necessary for the various types of calculations */
	switch (mgparm->type) {
		case MCT_MANUAL:
			return NOsh_setupCalcMGMANUAL(thee, calc);
		case MCT_DUMMY:
			return NOsh_setupCalcMGMANUAL(thee, calc);
		case MCT_AUTO:
			return NOsh_setupCalcMGAUTO(thee, calc);
		case MCT_PARALLEL:
			return NOsh_setupCalcMGPARA(thee, calc);
		default:
			Vnm_print(2, "NOsh_setupCalcMG:  undefined MG calculation type (%d)!\n",
					  mgparm->type);
			return 0;
	}
	
	/* Shouldn't get here */
	return 0;
}

VPRIVATE int NOsh_setupCalcFEM(
							   NOsh *thee,
							   NOsh_calc *calc
							   ) {
	
	VASSERT(thee != VNULL);
	VASSERT(calc != VNULL);
	VASSERT(calc->femparm != VNULL);
	
	/* Now we're ready to whatever sorts of post-processing operations that are
		* necessary for the various types of calculations */
	switch (calc->femparm->type) {
		case FCT_MANUAL:
			return NOsh_setupCalcFEMANUAL(thee, calc);
		default:
			Vnm_print(2, "NOsh_parseFEM:  unknown calculation type (%d)!\n",
					  calc->femparm->type);
			return 0;
	}
	
	/* Shouldn't get here */
	return 0;
}


VPRIVATE int NOsh_setupCalcMGMANUAL(
								   NOsh *thee, 
								   NOsh_calc *elec
								   ) {
	
	MGparm *mgparm = VNULL;
	PBEparm *pbeparm = VNULL;
	NOsh_calc *calc = VNULL;
	Valist *alist = VNULL;
	int i;
	
	if (thee == VNULL) {
		Vnm_print(2, "NOsh_setupCalcMGMANUAL:  Got NULL thee!\n");
		return 0;
	}
	if (elec == VNULL) {
		Vnm_print(2, "NOsh_setupCalcMGMANUAL:  Got NULL calc!\n");
		return 0;
	}
	mgparm = elec->mgparm;
	if (mgparm == VNULL) {
		Vnm_print(2, "NOsh_setupCalcMGMANUAL:  Got NULL mgparm -- was this calculation \
set up?\n");
		return 0;
	}
	pbeparm = elec->pbeparm;
	if (pbeparm == VNULL) {
		Vnm_print(2, "NOsh_setupCalcMGMANUAL:  Got NULL pbeparm -- was this calculation \
set up?\n");
		return 0;
	}
	
	/* Set up missing MG parameters */
	if (mgparm->setgrid == 0) {
		VASSERT(mgparm->setglen);
		mgparm->grid[0] = mgparm->glen[0]/((double)(mgparm->dime[0]-1));
		mgparm->grid[1] = mgparm->glen[1]/((double)(mgparm->dime[1]-1));
		mgparm->grid[2] = mgparm->glen[2]/((double)(mgparm->dime[2]-1));
	}
	if (mgparm->setglen == 0) {
		VASSERT(mgparm->setgrid);
		mgparm->glen[0] = mgparm->grid[0]*((double)(mgparm->dime[0]-1));
		mgparm->glen[1] = mgparm->grid[1]*((double)(mgparm->dime[1]-1));
		mgparm->glen[2] = mgparm->grid[2]*((double)(mgparm->dime[2]-1));	
	}
	
	/* Check to see if he have any room left for this type of calculation, if 
		so: set the calculation type, update the number of calculations of this type, 
		and parse the rest of the section */
	if (thee->ncalc >= NOSH_MAXCALC) {
		Vnm_print(2, "NOsh:  Too many calculations in this run!\n");
		Vnm_print(2, "NOsh:  Current max is %d; ignoring this calculation\n",
				  NOSH_MAXCALC);
		return 0;
	}
	
    /* Get the next calculation object and increment the number of calculations */
	thee->calc[thee->ncalc] = NOsh_calc_ctor(NCT_MG);
	calc = thee->calc[thee->ncalc];
	(thee->ncalc)++;
	
	/* Set up grid center on molecule if needed */
	if (mgparm->cmeth == MCM_MOLECULE) {
		VASSERT(mgparm->centmol >= 0);
		VASSERT(mgparm->centmol < thee->nmol);
		alist = thee->alist[mgparm->centmol];
		VASSERT(alist != VNULL);
		for (i=0; i<3; i++) {
			mgparm->center[i] = alist->center[i];
		}
	}
	
	/* Copy over contents of ELEC */
	NOsh_calc_copy(calc, elec);
	
	
    return 1;
}

VPUBLIC int NOsh_setupCalcMGAUTO(
								 NOsh *thee,
								 NOsh_calc *elec
								 ) {
	
	NOsh_calc *calcf = VNULL;
	NOsh_calc *calcc = VNULL;
	Valist *alist = VNULL;
	double fgrid[3], cgrid[3];
	double d[3], minf[3], maxf[3], minc[3], maxc[3];
	double redfrac, redrat[3], td, fshift, dshift;
	int ifocus, nfocus, tnfocus[3];
	int j;
	int icalc;
	int dofix;
	int molid;
	
	/* A comment about the coding style in this function.  I use lots and lots
		and lots of pointer deferencing.  I could (and probably should) save 
		these in temporary variables.  However, since there are so many MGparm, 
		etc. and NOsh_calc, etc. objects running around in this function, the 
		current scheme is easiest to debug. */

	
	if (thee == VNULL) {
		Vnm_print(2, "NOsh_setupCalcMGAUTO:  Got NULL thee!\n");
		return 0;
	}
	if (elec == VNULL) {
		Vnm_print(2, "NOsh_setupCalcMGAUTO:  Got NULL elec!\n");
		return 0;
	}
	if (elec->mgparm == VNULL) {
		Vnm_print(2, "NOsh_setupCalcMGAUTO:  Got NULL mgparm!\n");
		return 0;
	}
	if (elec->pbeparm == VNULL) {
		Vnm_print(2, "NOsh_setupCalcMGAUTO:  Got NULL pbeparm!\n");
		return 0;
	}
	
	/* Set elec coarse/fine centers to molecule, if requested */
	if (elec->mgparm->cmeth == MCM_MOLECULE) {
		molid = elec->mgparm->centmol;
		VASSERT(molid >= 0);
		VASSERT(molid < thee->nmol);
		alist = thee->alist[molid];
		VASSERT(alist != VNULL);
		for (j=0; j<3; j++) {
			elec->mgparm->center[j] = alist->center[j];
		}
	}
	if (elec->mgparm->fcmeth == MCM_MOLECULE) {
		molid = elec->mgparm->fcentmol;
		VASSERT(molid >= 0);
		VASSERT(molid < thee->nmol);
		alist = thee->alist[molid];
		VASSERT(alist != VNULL);
		for (j=0; j<3; j++) {
			if (elec->mgparm->type = MCT_PARALLEL) {
				elec->mgparm->fcenter[j] += alist->center[j];
			} else {
				elec->mgparm->fcenter[j] = alist->center[j];
			}
		}
	}
	if (elec->mgparm->ccmeth == MCM_MOLECULE) {
		molid = elec->mgparm->ccentmol;
		VASSERT(molid >= 0);
		VASSERT(molid < thee->nmol);
		alist = thee->alist[molid];
		VASSERT(alist != VNULL);
		for (j=0; j<3; j++) {
			elec->mgparm->ccenter[j] = alist->center[j];
		}
	}
	
	Vnm_print(0, "NOsh_setupCalcMGAUTO(%s, %d):  coarse grid center = %g %g %g\n",
			  __FILE__, __LINE__, 
			  elec->mgparm->ccenter[0],
			  elec->mgparm->ccenter[1],
			  elec->mgparm->ccenter[2]);
	Vnm_print(0, "NOsh_setupCalcMGAUTO(%s, %d):  fine grid center = %g %g %g\n",
			  __FILE__, __LINE__, 
			  elec->mgparm->fcenter[0],
			  elec->mgparm->fcenter[1],
			  elec->mgparm->fcenter[2]);
	
	/* Calculate the grid spacing on the coarse and fine levels */
	for (j=0; j<3; j++) {
		cgrid[j] = (elec->mgparm->cglen[j])/((double)(elec->mgparm->dime[j]-1));
		fgrid[j] = (elec->mgparm->fglen[j])/((double)(elec->mgparm->dime[j]-1));
		d[j] = elec->mgparm->fcenter[j] - elec->mgparm->ccenter[j];
	}
	Vnm_print(0, "NOsh:  Coarse grid spacing = %g, %g, %g\n", 
			  cgrid[0], cgrid[1], cgrid[2]);
	Vnm_print(0, "NOsh:  Fine grid spacing = %g, %g, %g\n",
			  fgrid[0], fgrid[1], fgrid[2]);
	Vnm_print(0, "NOsh:  Displacement between fine and coarse grids = %g, %g, %g\n",
			  d[0], d[1], d[2]);
	
	/* Now calculate the number of focusing levels, never reducing the grid 
		spacing by more than redfrac at each level */
	for (j=0; j<3; j++) {
		if (fgrid[j]/cgrid[j] < VREDFRAC) {
			redfrac = fgrid[j]/cgrid[j];
			td = log(redfrac)/log(VREDFRAC);
			tnfocus[j] = (int)ceil(td) + 1;
		} else tnfocus[j] = 2;
	}
	nfocus = VMAX2(VMAX2(tnfocus[0], tnfocus[1]), tnfocus[2]);
	
	/* Now set redrat to the actual value by which the grid spacing is reduced 
		at each level of focusing */
	for (j=0; j<3; j++) {
		redrat[j] = VPOW((fgrid[j]/cgrid[j]), 1.0/((double)nfocus-1.0));
	}
	Vnm_print(0, "NOsh:  %d levels of focusing with %g, %g, %g reductions\n",
			  nfocus, redrat[0], redrat[1], redrat[2]);
	
	/* Now that we know how many focusing levels to use, we're ready to set up 
		the parameter objects */
	if (nfocus > (NOSH_MAXCALC-(thee->ncalc))) {
		Vnm_print(2, "NOsh:  Require more calculations than max (%d)!\n",
				  NOSH_MAXCALC);
		return 0;
	}	

	for (ifocus=0; ifocus<nfocus; ifocus++) {
		
		/* Generate the new calc object */
		icalc = thee->ncalc;
		thee->calc[icalc] = NOsh_calc_ctor(NCT_MG);
		(thee->ncalc)++;
		
		/* This is the _current_ NOsh_calc object */
		calcf = thee->calc[icalc];
		/* This is the _previous_ Nosh_calc object */
		if (ifocus != 0) {
			calcc = thee->calc[icalc-1];
		} else {
			calcc = VNULL;
		}
		
		/* Copy over most of the parameters from the ELEC object */
		NOsh_calc_copy(calcf, elec);
		
		/* Set up the grid lengths and spacings */
		if (ifocus == 0) {
			for (j=0; j<3; j++) {
				calcf->mgparm->grid[j] = cgrid[j];
				calcf->mgparm->glen[j] = elec->mgparm->cglen[j];
			}
		} else {
			for (j=0; j<3; j++) {
				calcf->mgparm->grid[j] = redrat[j]*(calcc->mgparm->grid[j]);
                calcf->mgparm->glen[j] = redrat[j]*(calcc->mgparm->glen[j]);
            }
        }
        calcf->mgparm->setgrid = 1;
        calcf->mgparm->setglen = 1;
		
		/* Get centers and centering method from coarse and fine meshes */
        if (ifocus == 0) {
			calcf->mgparm->cmeth = elec->mgparm->ccmeth;
			calcf->mgparm->centmol = elec->mgparm->ccentmol;
			for (j=0; j<3; j++) {
				calcf->mgparm->center[j] = elec->mgparm->ccenter[j];
			}
		} else if (ifocus == (nfocus-1)) {
			calcf->mgparm->cmeth = elec->mgparm->fcmeth;
			calcf->mgparm->centmol = elec->mgparm->fcentmol;
			for (j=0; j<3; j++) {
				calcf->mgparm->center[j] = elec->mgparm->fcenter[j];
			}
		} else {
			calcf->mgparm->cmeth = MCM_FOCUS;
			/* TEMPORARILY move the current grid center 
			to the fine grid center.  In general, this will move portions of 
			the current mesh off the immediately-coarser mesh.  We'll fix that
			in the next step. */
			for (j=0; j<3; j++) {
				calcf->mgparm->center[j] = elec->mgparm->fcenter[j];
			}
		}

		
		/* As mentioned above, it is highly likely that the previous "jump" 
			to the fine grid center put portions of the current mesh off the 
			previous (coarser) mesh.  Fix this by displacing the current mesh 
			back onto the previous coarser mesh.  */
#define FUDGE_FRACTION 0.0
		if (ifocus != 0) {
			for (j=0; j<3; j++) {
				/* Check if we've fallen off of the lower end of the mesh */
				dofix = 0;
				minf[j] = calcf->mgparm->center[j] 
					- 0.5*(calcf->mgparm->glen[j]);
				minc[j] = calcc->mgparm->center[j] 
					- 0.5*(calcc->mgparm->glen[j]);
				d[j] = minc[j] - minf[j];
				if (d[j] >= VSMALL) {
					if (ifocus == (nfocus-1)) {
						Vnm_tprint(2, "NOsh_setupCalcMGAUTO:  Error!  Finest \
mesh has fallen off the coarser meshes!\n");
						Vnm_print(2, "NOsh_setupCalcMGAUTO:  difference in %d-\
direction = %g\n", j, d[j]);
						VASSERT(0);
					} else {
						Vnm_print(0, "NOsh_setupCalcMGAUTO(%s, %d):  ifocus = %d, \
fixing mesh min violation (%g in %d-direction).\n", __FILE__, __LINE__, ifocus, 
								  d[j], j);
						Vnm_print(0, "NOsh_setupCalcMGAUTO(%s, %d):  mesh center \
before fix = %g %g %g.\n", __FILE__, __LINE__, 
								  calcf->mgparm->center[0],
								  calcf->mgparm->center[1],
								  calcf->mgparm->center[2]
								  );
						calcf->mgparm->center[j] += \
						(d[j] + FUDGE_FRACTION * calcf->mgparm->grid[j]);
						Vnm_print(0, "NOsh_setupCalcMGAUTO(%s, %d):  mesh center \
after fix = %g %g %g.\n", __FILE__, __LINE__, 
								  calcf->mgparm->center[0],
								  calcf->mgparm->center[1],
								  calcf->mgparm->center[2]
								  );
						dofix = 1;
					}
				}
				/* Check if we've fallen off of the upper end of the mesh */
				maxf[j] = calcf->mgparm->center[j] \
					+ 0.5*(calcf->mgparm->glen[j]);
				maxc[j] = calcc->mgparm->center[j] \
					+ 0.5*(calcc->mgparm->glen[j]);
				d[j] = maxf[j] - maxc[j];
				if (d[j] >= VSMALL) {
					if (ifocus == (nfocus-1)) {
						Vnm_print(2, "NOsh_setupCalcMGAUTO:  Error!  Finest \
mesh has fallen off the coarser meshes!\n");
						Vnm_print(2, "NOsh_setupCalcMGAUTO:  difference in %d-\
direction = %g\n", j, d[j]);
						VASSERT(0);
					} else {
						/* If we already fixed the lower boundary and we now need
						to fix the upper boundary, we have a serious problem. */
						if (dofix) {
							Vnm_print(2, "NOsh_setupCalcMGAUTO:  Error!  Both \
ends of the finer mesh do not fit in the bigger mesh!\n");
							VASSERT(0);
						}
						Vnm_print(0, "NOsh_setupCalcMGAUTO(%s, %d):  ifocus = %d, \
fixing mesh max violation (%g in %d-direction).\n", __FILE__, __LINE__, ifocus, 
								  d[j], j);
						Vnm_print(0, "NOsh_setupCalcMGAUTO(%s, %d):  mesh center \
before fix = %g %g %g.\n", __FILE__, __LINE__, 
								  calcf->mgparm->center[0],
								  calcf->mgparm->center[1],
								  calcf->mgparm->center[2]
								  );
						calcf->mgparm->center[j] -= \
						(d[j] + FUDGE_FRACTION * calcf->mgparm->grid[j]);
						Vnm_print(0, "NOsh_setupCalcMGAUTO(%s, %d):  mesh center \
after fix = %g %g %g.\n", __FILE__, __LINE__, 
								  calcf->mgparm->center[0],
								  calcf->mgparm->center[1],
								  calcf->mgparm->center[2]
								  );
						dofix = 1;
					}
				}
			}
		}	
#undef FUDGE_FRACTION
		
		/* Finer levels have focusing boundary conditions */
		if (ifocus != 0) calcf->pbeparm->bcfl = BCFL_FOCUS;
		
		/* Only the finest level handles I/O and needs to worry about disjoint
			partitioning */
		if (ifocus != (nfocus-1)) calcf->pbeparm->numwrite = 0;
		
		/* Reset boundary flags for parallel focusing */ 
        if (calcf->mgparm->type != 2)  {
			Vnm_print(0, "NOsh_setupMGAUTO:  Resetting boundary flags\n");
			for (j=0; j<6; j++) calcf->mgparm->partDisjOwnSide[j] = 0;
			for (j=0; j<3; j++) {
				calcf->mgparm->partDisjCenterShift[j] = 0;
				calcf->mgparm->partDisjLength[j] = calcf->mgparm->glen[j];
			}
        }
		
		
        calcf->mgparm->parsed = 1;
    }
	
	
    return 1;
}

/* Author:   Nathan Baker and Todd Dolinsky */
VPUBLIC int NOsh_setupCalcMGPARA(
								 NOsh *thee, 
								 NOsh_calc *elec
								 ) {
	
	/* NEW (25-Jul-2006):  This code should produce modify the ELEC statement
	and pass it on to MGAUTO for further processing. */
	
	MGparm *mgparm = VNULL;
	double ofrac;
	double hx, hy, hzed;
	double xofrac, yofrac, zofrac;
	int rank, size, npx, npy, npz, nproc, ip, jp, kp;
	int xeffGlob, yeffGlob, zeffGlob, xDisj, yDisj, zDisj;
	int xigminDisj, xigmaxDisj, yigminDisj, yigmaxDisj, zigminDisj, zigmaxDisj;
	int xigminOlap, xigmaxOlap, yigminOlap, yigmaxOlap, zigminOlap, zigmaxOlap;
	int xOlapReg, yOlapReg, zOlapReg;
	double xlenDisj, ylenDisj, zlenDisj;
	double xcentDisj, ycentDisj, zcentDisj;
	double xcentOlap, ycentOlap, zcentOlap;
	double xlenOlap, ylenOlap, zlenOlap;
	double xminOlap, xmaxOlap, yminOlap, ymaxOlap, zminOlap, zmaxOlap;
	double xminDisj, xmaxDisj, yminDisj, ymaxDisj, zminDisj, zmaxDisj;
	double xcent, ycent, zcent;
	
	/* Grab some useful variables */
	VASSERT(thee != VNULL);
	VASSERT(elec != VNULL);
	mgparm = elec->mgparm;
	VASSERT(mgparm != VNULL);
	
	/* Grab some useful variables */
    ofrac = mgparm->ofrac;  
    npx = mgparm->pdime[0]; 
    npy = mgparm->pdime[1]; 
    npz = mgparm->pdime[2]; 
    nproc = npx*npy*npz;
	
	/* If this is not an asynchronous calculation, then we need to make sure we
		have all the necessary MPI information */
    if (mgparm->setasync == 0) {
		
#ifndef HAVE_MPI_H
		
        Vnm_tprint(2, "NOsh_setupCalcMGPARA:  Oops!  You're trying to perform \
an 'mg-para' (parallel) calculation\n");
        Vnm_tprint(2, "NOsh_setupCalcMGPARA:  with a version of APBS that wasn't \
compiled with MPI!\n");
		Vnm_tprint(2, "NOsh_setupCalcMGPARA:  Perhaps you meant to use the \
'async' flag?\n");
        Vnm_tprint(2, "NOsh_setupCalcMGPARA:  Bailing out!\n");
		
        return 0;
		
#endif
		
		rank = thee->proc_rank;
		size = thee->proc_size;
		Vnm_print(0, "NOsh_setupCalcMGPARA:  Hello from processor %d of %d\n", rank,
				  size);
		
		/* Check to see if we have too many processors.  If so, then simply set
			this processor to duplicating the work of processor 0. */
		if (rank > (nproc-1)) {
			Vnm_print(2, "NOsh_setupMGPARA:  There are more processors available than\
the %d you requested.\n", nproc);
			Vnm_print(2, "NOsh_setupMGPARA:  Eliminating processor %d\n", rank);
			thee->bogus = 1;
			rank = 0;
		}
		
		/* Check to see if we have too few processors.  If so, this is a fatal 
			error. */
		if (size < nproc) {
			Vnm_print(2, "NOsh_setupMGPARA:  There are too few processors (%d) to \
satisfy requirements (%d)\n", size, nproc);
			return 0;
		}
		
		Vnm_print(0, "NOsh_setupMGPARA:  Hello (again) from processor %d of %d\n", 
				  rank, size);
		
    } else { /* Setting up for an asynchronous calculation. */
		
        rank = mgparm->async;
		
        /* Check to see if the async id is greater than the number of
		* processors.  If so, this is a fatal error. */
        if (rank > (nproc-1)) {
            Vnm_print(2, "NOsh_setupMGPARA:  The processor id you requested (%d) \
is not within the range of processors available (0-%d)\n", rank, (nproc-1));
            return 0;
        }
    }
	

	
    /* Calculate the processor's coordinates in the processor grid */
	kp = (int)floor(rank/(npx*npy));
	jp = (int)floor((rank-kp*npx*npy)/npx);
	ip = rank - kp*npx*npy - jp*npx;
	Vnm_print(0, "NOsh_setupMGPARA:  Hello world from PE (%d, %d, %d)\n", 
			  ip, jp, kp);
	
	/* Calculate effective overlap fractions for uneven processor distributions */
	if (npx == 1) xofrac = 0.0;
	else xofrac = ofrac;
	if (npy == 1) yofrac = 0.0;
	else yofrac = ofrac;
	if (npz == 1) zofrac = 0.0;
	else zofrac = ofrac;
	
    /* Calculate the global grid size and spacing */
    xDisj = (int)VFLOOR(mgparm->dime[0]/(1 + 2*xofrac) + 0.5);
    xeffGlob = npx*xDisj;
    hx = mgparm->fglen[0]/(double)(xeffGlob-1);
    yDisj = (int)VFLOOR(mgparm->dime[1]/(1 + 2*yofrac) + 0.5);
    yeffGlob = npy*yDisj;
    hy = mgparm->fglen[1]/(double)(yeffGlob-1);
    zDisj = (int)VFLOOR(mgparm->dime[2]/(1 + 2*zofrac) + 0.5);
    zeffGlob = npz*zDisj;
    hzed = mgparm->fglen[2]/(double)(zeffGlob-1);
    Vnm_print(0, "NOsh_setupMGPARA:  Global Grid size = (%d, %d, %d)\n",
              xeffGlob, yeffGlob, zeffGlob);
    Vnm_print(0, "NOsh_setupMGPARA:  Global Grid Spacing = (%.3f, %.3f, %.3f)\n",
              hx, hy, hzed);
    Vnm_print(0, "NOsh_setupMGPARA:  Processor Grid Size = (%d, %d, %d)\n",
              xDisj, yDisj, zDisj);
	
    /* Calculate the maximum and minimum processor grid points */	
    xigminDisj = ip*xDisj;
    xigmaxDisj = xigminDisj + xDisj - 1;
    yigminDisj = jp*yDisj;
    yigmaxDisj = yigminDisj + yDisj - 1;
    zigminDisj = kp*zDisj;
    zigmaxDisj = zigminDisj + zDisj - 1;
    Vnm_print(0, "NOsh_setupMGPARA:  Min Grid Points for this proc. (%d, %d, %d)\n",
              xigminDisj, yigminDisj, zigminDisj);
    Vnm_print(0, "NOsh_setupMGPARA:  Max Grid Points for this proc. (%d, %d, %d)\n",
              xigmaxDisj, yigmaxDisj, zigmaxDisj);
	
	
    /* Calculate the disjoint partition length and center displacement */
    xminDisj = VMAX2(hx*(xigminDisj-0.5), 0.0);
    xmaxDisj = VMIN2(hx*(xigmaxDisj+0.5), mgparm->fglen[0]);
    xlenDisj = xmaxDisj - xminDisj;
    yminDisj = VMAX2(hy*(yigminDisj-0.5), 0.0);
    ymaxDisj = VMIN2(hy*(yigmaxDisj+0.5), mgparm->fglen[1]);
    ylenDisj = ymaxDisj - yminDisj;
    zminDisj = VMAX2(hzed*(zigminDisj-0.5), 0.0);
    zmaxDisj = VMIN2(hzed*(zigmaxDisj+0.5), mgparm->fglen[2]);
    zlenDisj = zmaxDisj - zminDisj;       
	
    xcent = 0.5*mgparm->fglen[0];
    ycent = 0.5*mgparm->fglen[1];
    zcent = 0.5*mgparm->fglen[2];
    
    xcentDisj = xminDisj + 0.5*xlenDisj - xcent;
    ycentDisj = yminDisj + 0.5*ylenDisj - ycent;
    zcentDisj = zminDisj + 0.5*zlenDisj - zcent;
    if (VABS(xcentDisj) < VSMALL) xcentDisj = 0.0;
    if (VABS(ycentDisj) < VSMALL) ycentDisj = 0.0;
    if (VABS(zcentDisj) < VSMALL) zcentDisj = 0.0;
	
    Vnm_print(0, "NOsh_setupMGPARA:  Disj part length = (%g, %g, %g)\n", 
			  xlenDisj, ylenDisj, zlenDisj);
    Vnm_print(0, "NOsh_setupMGPARA:  Disj part center displacement = (%g, %g, %g)\n", 
			  xcentDisj, ycentDisj, zcentDisj);
	
    /* Calculate the overlapping partition length and center displacement */
    xOlapReg = 0;
    yOlapReg = 0;
    zOlapReg = 0;
    if (npx != 1) xOlapReg = (int)VFLOOR(xofrac*mgparm->fglen[0]/npx/hx + 0.5) + 1;
    if (npy != 1) yOlapReg = (int)VFLOOR(yofrac*mgparm->fglen[1]/npy/hy + 0.5) + 1;
    if (npz != 1) zOlapReg = (int)VFLOOR(zofrac*mgparm->fglen[2]/npz/hzed + 0.5) + 1;
	
    Vnm_print(0, "NOsh_setupMGPARA:  No. of Grid Points in Overlap (%d, %d, %d)\n",
              xOlapReg, yOlapReg, zOlapReg);
	
    if (ip == 0) xigminOlap = 0;
    else if (ip == (npx - 1)) xigminOlap = xeffGlob - mgparm->dime[0];
    else xigminOlap = xigminDisj - xOlapReg;
    xigmaxOlap = xigminOlap + mgparm->dime[0] - 1;
	
    if (jp == 0) yigminOlap = 0;
    else if (jp == (npy - 1)) yigminOlap = yeffGlob - mgparm->dime[1];
    else yigminOlap = yigminDisj - yOlapReg;
    yigmaxOlap = yigminOlap + mgparm->dime[1] - 1;
	
    if (kp == 0) zigminOlap = 0;
    else if (kp == (npz - 1)) zigminOlap = zeffGlob - mgparm->dime[2];
    else zigminOlap = zigminDisj - zOlapReg;
    zigmaxOlap = zigminOlap + mgparm->dime[2] - 1;
    
    Vnm_print(0, "NOsh_setupMGPARA:  Min Grid Points with Overlap (%d, %d, %d)\n",
              xigminOlap, yigminOlap, zigminOlap);
    Vnm_print(0, "NOsh_setupMGPARA:  Max Grid Points with Overlap (%d, %d, %d)\n",
              xigmaxOlap, yigmaxOlap, zigmaxOlap);
	
    xminOlap = hx * xigminOlap;
    xmaxOlap = hx * xigmaxOlap;
    yminOlap = hy * yigminOlap;
    ymaxOlap = hy * yigmaxOlap;
    zminOlap = hzed * zigminOlap;
    zmaxOlap = hzed * zigmaxOlap;
	
    xlenOlap = xmaxOlap - xminOlap;
    ylenOlap = ymaxOlap - yminOlap;
    zlenOlap = zmaxOlap - zminOlap;
	
    xcentOlap = (xminOlap + 0.5*xlenOlap) - xcent;
    ycentOlap = (yminOlap + 0.5*ylenOlap) - ycent;
    zcentOlap = (zminOlap + 0.5*zlenOlap) - zcent;
    if (VABS(xcentOlap) < VSMALL) xcentOlap = 0.0;
    if (VABS(ycentOlap) < VSMALL) ycentOlap = 0.0;
    if (VABS(zcentOlap) < VSMALL) zcentOlap = 0.0;
    
    Vnm_print(0, "NOsh_setupMGPARA:  Olap part length = (%g, %g, %g)\n", 
			  xlenOlap, ylenOlap, zlenOlap);
    Vnm_print(0, "NOsh_setupMGPARA:  Olap part center displacement = (%g, %g, %g)\n", 
			  xcentOlap, ycentOlap, zcentOlap);
	
	
    /* Calculate the boundary flags:
		Flags are set to 1 when another processor is present along the boundary
		Flags are otherwise set to 0. */
	
    if (ip == 0) mgparm->partDisjOwnSide[VAPBS_LEFT] = 0;
	else mgparm->partDisjOwnSide[VAPBS_LEFT] = 1;
	if (ip == (npx-1)) mgparm->partDisjOwnSide[VAPBS_RIGHT] = 0;
	else mgparm->partDisjOwnSide[VAPBS_RIGHT] = 1;
	if (jp == 0) mgparm->partDisjOwnSide[VAPBS_BACK] = 0;
	else mgparm->partDisjOwnSide[VAPBS_BACK] = 1;
	if (jp == (npy-1)) mgparm->partDisjOwnSide[VAPBS_FRONT] = 0;
	else mgparm->partDisjOwnSide[VAPBS_FRONT] = 1;
	if (kp == 0) mgparm->partDisjOwnSide[VAPBS_DOWN] = 0;
	else mgparm->partDisjOwnSide[VAPBS_DOWN] = 1;
	if (kp == (npz-1)) mgparm->partDisjOwnSide[VAPBS_UP] = 0;
	else mgparm->partDisjOwnSide[VAPBS_UP] = 1;       
	
	Vnm_print(0, "NOsh_setupMGPARA:  partDisjOwnSide[LEFT] = %d\n", 
			  mgparm->partDisjOwnSide[VAPBS_LEFT]);
	Vnm_print(0, "NOsh_setupMGPARA:  partDisjOwnSide[RIGHT] = %d\n", 
			  mgparm->partDisjOwnSide[VAPBS_RIGHT]);
	Vnm_print(0, "NOsh_setupMGPARA:  partDisjOwnSide[FRONT] = %d\n", 
			  mgparm->partDisjOwnSide[VAPBS_FRONT]);
	Vnm_print(0, "NOsh_setupMGPARA:  partDisjOwnSide[BACK] = %d\n", 
			  mgparm->partDisjOwnSide[VAPBS_BACK]);
	Vnm_print(0, "NOsh_setupMGPARA:  partDisjOwnSide[UP] = %d\n", 
			  mgparm->partDisjOwnSide[VAPBS_UP]);
	Vnm_print(0, "NOsh_setupMGPARA:  partDisjOwnSide[DOWN] = %d\n", 
			  mgparm->partDisjOwnSide[VAPBS_DOWN]);
	
	
	/* Set the disjoint partition information */
	mgparm->partDisjCenterShift[0] = xcentDisj;
	mgparm->partDisjCenterShift[1] = ycentDisj;
	mgparm->partDisjCenterShift[2] = zcentDisj;
	mgparm->partDisjLength[0] = xlenDisj;
	mgparm->partDisjLength[1] = ylenDisj;
	mgparm->partDisjLength[2] = zlenDisj;
	
	/* Set the fine mesh parameters to the non-disjoint partition */ 
	mgparm->fglen[0] = xlenOlap;
	mgparm->fglen[1] = ylenOlap;
	mgparm->fglen[2] = zlenOlap;
	mgparm->fcenter[0] += xcentOlap;
	mgparm->fcenter[1] += ycentOlap;
	mgparm->fcenter[2] += zcentOlap;
	
	/* Setup the automatic focusing calculations associated with this processor */
	return NOsh_setupCalcMGAUTO(thee, elec);

}

VPUBLIC int NOsh_parseFEM(
						  NOsh *thee, 
						  Vio *sock, 
						  NOsh_calc *elec
						  ) {
	
	char tok[VMAX_BUFSIZE];
	FEMparm *feparm = VNULL;    
	PBEparm *pbeparm = VNULL;    
	int rc;
	
	/* Check the arguments */
	if (thee == VNULL) {
		Vnm_print(2, "NOsh_parseFEM:  Got NULL thee!\n");
		return 0;
	}
	if (sock == VNULL) {
		Vnm_print(2, "NOsh_parseFEM:  Got pointer to NULL socket!\n");
		return 0;
	}
	if (elec == VNULL) {
		Vnm_print(2, "NOsh_parseFEM:  Got pointer to NULL elec object!\n");
		return 0;
	}
	feparm = elec->femparm;
	if (feparm == VNULL) {
		Vnm_print(2, "NOsh_parseFEM:  Got pointer to NULL feparm object!\n");
		return 0;
	}
	pbeparm = elec->pbeparm;
	if (feparm == VNULL) {
		Vnm_print(2, "NOsh_parseFEM:  Got pointer to NULL pbeparm object!\n");
		return 0;
	}
	
	Vnm_print(0, "NOsh_parseFEM: Parsing parameters for FEM calculation\n");

	/* Start snarfing tokens from the input stream */
	rc = 1;
	while (Vio_scanf(sock, "%s", tok) == 1) { 
		
		Vnm_print(0, "NOsh_parseFEM:  Parsing %s...\n", tok);
		
		/* See if it's an END token */
		if (Vstring_strcasecmp(tok, "end") == 0) {
			feparm->parsed = 1;
			pbeparm->parsed = 1;
			rc = 1;
			break;
		}
		
		/* Pass the token through a series of parsers */
		rc = PBEparm_parseToken(pbeparm, tok, sock);
		if (rc == -1) {
			Vnm_print(0, "NOsh_parseFEM:  parsePBE error!\n");
			break;
		} else if (rc == 0) {
			/* Pass the token to the generic MG parser */
			rc = FEMparm_parseToken(feparm, tok, sock);
			if (rc == -1) { 
				Vnm_print(0, "NOsh_parseFEM:  parseMG error!\n");
				break;
			} else if (rc == 0) {
				/* We ran out of parsers! */
				Vnm_print(2, "NOsh:  Unrecognized keyword: %s\n", tok);
				break;
			}
		}
	}
	
	/* Handle various errors arising in the token-snarfing loop -- these all
		* just result in simple returns right now */
	if (rc == -1) return 0;
	if (rc == 0) return 0;
	
	/* Check the status of the parameter objects */
	if ((!FEMparm_check(feparm)) || (!PBEparm_check(pbeparm))) {
		Vnm_print(2, "NOsh:  FEM parameters not set correctly!\n");
		return 0;
	}
	
	return 1;
	
}

VPRIVATE int NOsh_setupCalcFEMANUAL(
								   NOsh *thee, 
								   NOsh_calc *elec
								   ) {
	
	FEMparm *feparm = VNULL;
	PBEparm *pbeparm = VNULL;
	NOsh_calc *calc = VNULL;
	Valist *alist = VNULL;
	int i;
	
	VASSERT(thee != VNULL);
	VASSERT(elec != VNULL);
	feparm = elec->femparm;
	VASSERT(feparm != VNULL);
	pbeparm = elec->pbeparm;
	VASSERT(pbeparm);
	
	/* Check to see if he have any room left for this type of
	* calculation, if so: set the calculation type, update the number
	* of calculations of this type, and parse the rest of the section
	*/
	if (thee->ncalc >= NOSH_MAXCALC) {
		Vnm_print(2, "NOsh:  Too many calculations in this run!\n");
		Vnm_print(2, "NOsh:  Current max is %d; ignoring this calculation\n",
				  NOSH_MAXCALC);
		return 0;
	}
	thee->calc[thee->ncalc] = NOsh_calc_ctor(NCT_FEM);
	calc = thee->calc[thee->ncalc];
	(thee->ncalc)++;
	
	/* Copy over contents of ELEC */
	NOsh_calc_copy(calc, elec);
	
	
	return 1;
}
