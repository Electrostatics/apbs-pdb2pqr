/**
 *  @file    benchmark.c
 *  @author  Nathan Baker
 *  @brief   Sample program for benchmarking I/O
 *  @version $Id$
 *
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
 * Copyright (c) 2002-2005.  Washington University in St. Louis.
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
 * code included in releases of ISIM, PMV, PyMOL, SMOL, VMD, and Vision.
 * Such combined software may be linked with APBS and redistributed together 
 * in original or modified form as mere aggregation without requirement that 
 * the entire work be under the scope of the GNU General Public License.
 * This special exception permission is also extended to any software listed
 * in the SPECIAL GPL EXCEPTION clauses by the PMG, FEtk, MC, or MALOC
 * libraries.
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

#define ERRRC 12

VEMBED(rcsid="$Id$")

int main(int argc, char **argv) {

    /* *************** VARIABLES ******************* */
    int nx, ny, nz, i, j, k, itmp;
    double hy, hx, hzed, xmin, ymin, zmin, pt[3], value, sum, dtmp;
    int test_int = 44;
    double test_double = -37.56;
    double test_exp = 5555000.12;
    char *usage = "\n  benchmark <asc|xdr> nx ny nz\n"; 
    char *iofmt;
    char test_string[VMAX_BUFSIZE];
    char tok[VMAX_BUFSIZE];
    char *MCwhiteChars = " =,;\t\n";
    char *MCcommChars  = "#%";
    Vgrid *grid;
    Vio *sock;
 
    /* *************** CHECK INVOCATION ******************* */
    Vio_start();
    Vnm_redirect(0);
    if (argc != 5) {
        Vnm_print(2,"\n*** Syntax error: got %d arguments, expected 4.\n\n",
          argc);
        Vnm_print(2,"%s\n", usage);
        return ERRRC;
    }
    if (Vstring_strcasecmp(argv[1], "XDR") == 0) iofmt = "XDR";
    else if (Vstring_strcasecmp(argv[1], "ASC") == 0) iofmt = "ASC";
    else { 
        Vnm_print(2, "Invalid format (%s)!\n", argv[1]);
        Vnm_print(2,"%s\n", usage);
        return ERRRC;
    }
    sscanf(argv[2], "%d", &nx);
    sscanf(argv[3], "%d", &ny);
    sscanf(argv[4], "%d", &nz);
    hx = 1;  
    hy = 1;
    hzed = 1;
    xmin = 0;  
    ymin = 0;
    zmin = 0;

    Vnm_print(1, "Validating I/O for %s format...\n", iofmt);
    sock = Vio_ctor("FILE", iofmt, "localhost", "benchmark.test", "w");
    if (sock == VNULL) {
        Vnm_print(2, "Problem opening virtual socket!\n");
        return ERRRC;
    }
    if (Vio_connect(sock, 0) < 0) {
        Vnm_print(2, "Vgrid_readDX: Problem connecting virtual socket!\n");
        return ERRRC;
    }
    Vio_setWhiteChars(sock, MCwhiteChars);
    Vio_setCommChars(sock, MCcommChars);
    sprintf(test_string, "integer %d double %4.3f exponential %12.5e\n", 
      test_int, test_double, test_exp);
    Vnm_print(1, "Writing '%s' to socket\n", test_string);
    Vio_printf(sock, "%s", test_string);
    Vio_connectFree(sock);
    Vio_dtor(&sock);
    
    sock = Vio_ctor("FILE", iofmt, "localhost", "benchmark.test", "r");
    if (sock == VNULL) {
        Vnm_print(2, "Problem opening virtual socket!\n");
        return ERRRC;
    }
    if (Vio_accept(sock, 0) < 0) {
        Vnm_print(2, "Vgrid_readDX: Problem accepting virtual socket!\n");
        return ERRRC;
    }
    Vio_setWhiteChars(sock, MCwhiteChars);
    Vio_setCommChars(sock, MCcommChars);
    VASSERT(Vio_scanf(sock, "%s", tok) == 1);
    Vnm_print(1, "Read token '%s' from socket\n", tok);
    VASSERT(Vio_scanf(sock, "%s", tok) == 1);
    VASSERT(sscanf(tok, "%i", &itmp) == 1);
    Vnm_print(1, "Read integer '%d' from socket\n", itmp);
    VASSERT(itmp == test_int);
    VASSERT(Vio_scanf(sock, "%s", tok) == 1);
    Vnm_print(1, "Read token '%s' from socket\n", tok);
    VASSERT(Vio_scanf(sock, "%s", tok) == 1);
    VASSERT(sscanf(tok, "%lf", &dtmp) == 1);
    Vnm_print(1, "Read double '%lf' from socket\n", dtmp);
    VASSERT(dtmp == test_double);
    VASSERT(Vio_scanf(sock, "%s", tok) == 1);
    Vnm_print(1, "Read token '%s' from socket\n", tok);
    VASSERT(Vio_scanf(sock, "%s", tok) == 1);
    VASSERT(sscanf(tok, "%le", &dtmp) == 1);
    Vnm_print(1, "Read exponential '%le' from socket\n", dtmp);
    VASSERT(dtmp == test_exp);
    Vio_acceptFree(sock);
    Vio_dtor(&sock);

    return 0;

}
