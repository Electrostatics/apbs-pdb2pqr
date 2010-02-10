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
 * Copyright (c) 2002-2010, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2010.  Nathan A. Baker
 * Portions Copyright (c) 1999-2002.  The Regents of the University of California.
 * Portions Copyright (c) 1995.  Michael Holst
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  
 * 
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. *
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
