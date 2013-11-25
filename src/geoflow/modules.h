///  @file    modules.h
///  @author  Andrew Stevens, Kyle Monson
///  @brief FORTRAN globals and function headers
///  @ingroup Geoflow
///  @version $Id$
///  @attention
///  @verbatim
///
/// APBS -- Adaptive Poisson-Boltzmann Solver
///
///  Nathan A. Baker (nathan.baker@pnnl.gov)
///  Pacific Northwest National Laboratory
///
///  Additional contributing authors listed in the code documentation.
///
/// Copyright (c) 2010-2012 Battelle Memorial Institute. Developed at the
/// Pacific Northwest National Laboratory, operated by Battelle Memorial
/// Institute, Pacific Northwest Division for the U.S. Department of Energy.
///
/// Portions Copyright (c) 2002-2010, Washington University in St. Louis.
/// Portions Copyright (c) 2002-2010, Nathan A. Baker.
/// Portions Copyright (c) 1999-2002, The Regents of the University of
/// California.
/// Portions Copyright (c) 1995, Michael Holst.
/// All rights reserved.
///
/// Redistribution and use in source and binary forms, with or without
/// modification, are permitted provided that the following conditions are met:
///
/// Redistributions of source code must retain the above copyright notice, this
/// list of conditions and the following disclaimer.
///
/// Redistributions in binary form must reproduce the above copyright notice,
/// this list of conditions and the following disclaimer in the documentation
/// and/or other materials provided with the distribution.
///
/// Neither the name of the developer nor the names of its contributors may be
/// used to endorse or promote products derived from this software without
/// specific prior written permission.
///
/// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
/// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
/// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
/// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
/// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
/// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
/// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
/// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
/// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
/// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
/// THE POSSIBILITY OF SUCH DAMAGE.
///
/// @endverbatim

#ifndef MODULES_H
#define MODULES_H

#include <vector>
#include "Mat.h"

#define MAXATOMS 15000  /* From the original f90 code; need to keep this
             * constant since the f90 routines that we call from here 
             * depend on it. */
#define XYZRWIDTH 4     // 4: 0-2 => pos, 3 => radius

#define MAX(x,y)  ((x > y) ? x : y)

struct Comdata{
    char fname[100];
    int nx, ny, nz;
    double xleft, xright,
           yleft, yright,
           zleft, zright,
    
           deltax, deltay, deltaz,
    
           dcel,
           pi;
    std::vector<double> xc, yc, zc;
};
extern Comdata comdata;

struct LJ{
    double tauval, prob, vdwdispersion,
           sigmas, roro, conms,
           density, epsilonw;
    int ffmodel;
    static const int iosetar = 1, iosetaa = 1, iwca = 1;
};
extern LJ lj;

double dot(double x, double y, double z);

extern "C"{
void domainini(double xyzr[MAXATOMS][XYZRWIDTH], const size_t natm, const double extvalue);

void chargedist(double xyzr[MAXATOMS][XYZRWIDTH], double* chratm, Mat<>& charget, Mat<>& corlocqt, Mat<int>& loc_qt, size_t iatm);

void yhsurface(double xyzr[MAXATOMS][XYZRWIDTH], double* ljepsilon, size_t natm, double tott,
    double deltat, Mat<>& phix, Mat<>& surfu, int i, double& area, double& vol, double& attint,
    double alpha, int iadi, int igfin);

void seteqb(Mat<>& bg, double xyzr[MAXATOMS][XYZRWIDTH], double* pqr, Mat<>& charget, Mat<>& corlocqt, double *epsilonsp);
//  void writerms_gama( double* sumpot, double* expv, double* elec, size_t* natm, double* gama, int *ngiter, double *pres );

void pbsolver( Mat<>& eps, Mat<>& phi, Mat<>& bg, double dcel,  double tol, int iter);

double xvalue(int& i);
double yvalue(int& i);
double zvalue(int& i);

int inverx(double& x);
int invery(double& y);
int inverz(double& z);
}

#endif

