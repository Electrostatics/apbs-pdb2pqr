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
/// PARTICULAR PURPOSE.  THE SOFTWARE AND PMGOMPANYING DOCUMENTATION, IF
/// ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
/// TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
/// MODIFICATIONS. 
////////////////////////////////////////////////////////////////////////////
/// rcsid="$Id$"
//////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// File:     vpmg-private.h    < vpmg.c >
//
// Purpose:  Private declarations for class Vpmg
//     
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VPMG_PRIVATE_H_
#define _VPMG_PRIVATE_H_

#include "apbscfg.h"

/* ///////////////////////////////////////////////////////////////////////////
// External FORTRAN ROUTINES 
/////////////////////////////////////////////////////////////////////////// */
#define F77MGSZ       VF77_MANGLE(mgsz, MGSZ)
#define F77PACKMG     VF77_MANGLE(packmg, PACKMG)
#define F77CGMGDRIV   VF77_MANGLE(cgmgdriv, CGMGDRIV)
#define F77NEWDRIV    VF77_MANGLE(newdriv, NEWDRIV)
#define F77MGDRIV     VF77_MANGLE(mgdriv, MGDRIV)
#define F77NCGHSDRIV  VF77_MANGLE(ncghsdriv, NCGHSDRIV)
#define F77NSORDRIV   VF77_MANGLE(nsordriv, NSORDRIV)
#define F77NGSRBDRIV  VF77_MANGLE(ngsrbdriv, NGSRBDRIV)
#define F77NWJACDRIV  VF77_MANGLE(nwjacdriv, NWJACDRIV)
#define F77NRICHDRIV  VF77_MANGLE(nrichdriv, NRICHDRIV)
#define F77TSECND     VF77_MANGLE(tsecnd, TSECND)
#define F77VPMGANORM  VF77_MANGLE(vpmganorm, VPMGANORM)
#define F77VPMGABAND  VF77_MANGLE(vpmgaband, VPMGABAND)
#define F77DPBFA      VF77_MANGLE(dpbfa, DPBFA)
#define F77DPBDI      VF77_MANGLE(dpbdi, DPBDI)
#define F77EIGDRIV    VF77_MANGLE(eigdriv, EIGDRIV)
#define F77ANORMDRIV  VF77_MANGLE(anormdriv, ANORMDRIV)
#define F77MYPDEFINIT  VF77_MANGLE(mypdefinit, MYPDEFINIT)

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmg: Private methods
/////////////////////////////////////////////////////////////////////////// */
#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define IJKx(j,k,i) (((i)*(ny)*(nz))+((k)*(ny))+(j))
#define IJKy(i,k,j) (((j)*(nx)*(nz))+((k)*(nx))+(i))
#define IJKz(i,j,k) (((k)*(nx)*(ny))+((j)*(nx))+(i))

#endif
