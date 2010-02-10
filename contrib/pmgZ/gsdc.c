/**
 *  @file    gsd.c
 *  @ingroup PMG
 *  @author  David Gohara
 *  @brief   C Extensions for FORTRAN PMG
 *  @version 
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
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MIN(a,b) ((a > b) ? b : a)
#define MAX(a,b) ((a > b) ? a : b)

void doKernel(int i1,int i2,int j1,int j2,int k1,int k2,int nx,int ny,int nz,
			  double * oN,double * oE,double * uC, 
			  double * oC,double * x,double * fc,double * cc,int red){
	
	int i,j,k;
	
	int kxy, jx, ioff;
	int ijk, ijp1k, ijm1k;
	int ip1jk, im1jk, ijkp1, ijkm1;
	
	double inter1, inter2;
	
	for(k=k1;k<k2;k++){
		kxy = k*nx*ny;
		for(j=j1;j<j2;j++){
			jx = j*nx;
			ioff = red ? ((j+k) % 2) : (1 - ((j+k) % 2));
			for(i=i1+ioff;i<i2;i+=2){
				ijk = kxy+jx+i;
				
				ijp1k = kxy+(j+1)*nx+i;
				ijm1k = kxy+(j-1)*nx+i;
				
				ip1jk = (k+1)*nx*ny+jx+i;
				im1jk = (k-1)*nx*ny+jx+i;
				
				ijkp1 = kxy+jx+(i+1);
				ijkm1 = kxy+jx+(i-1);
				
				inter1 = (fc[ijk] + (
									 + oN[ijk] * x[ijp1k]
									 + oN[ijm1k] * x[ijm1k]
									 + oE[ijk] * x[ijkp1]
									 + oE[ijkm1] * x[ijkm1]
									 + uC[ijk] * x[ip1jk]
									 + uC[im1jk] * x[im1jk]
									 ));
				inter2 = (oC[ijk] + cc[ijk]);
				
				x[ijk] =  inter1 / inter2;
			}
		}
	}
}

void gsrb7c_(double * oN,double * oE,double * uC, double * oC,
			 double * x,double * fc,double * cc,int * iadjoint,
			 int * aNx,int * aNy,int * aNz,int * aItmax,
			 double * r,int doResid){
	
	int i,j,k, iters;
	int i1,i2,j1,j2,k1,k2;
	int nx, ny, nz, itmax;
	int ti, tj, tk;
	
	itmax = *aItmax;
	
	nx = *aNx;
	ny = *aNy;
	nz = *aNz;
	
	ti = nx;
	tj = ny;
	tk = 1;
	
	for (iters=0;iters<itmax;iters++){	
//#pragma omp parallel for private(i,j,k,i1,i2,j1,j2,k1,k2) schedule(dynamic,tk)
		for(k=1;k<nz-1;k+=tk){
			k1 = k;
			k2 = MIN(k+tk,nz-1);
			for(j=1;j<ny-1;j+=tj){
				j1 = j;
				j2 = MIN(j+tj,ny-1);
				for(i=1;i<nx-1;i+=ti){
					i1 = i;
					i2 = MIN(i+ti,nx-1);
					doKernel(i1,i2,j1,j2,k1,k2,nx,ny,nz,oN,oE,uC,oC,x,fc,cc,1);
					doKernel(i1,i2,j1,j2,k1,k2,nx,ny,nz,oN,oE,uC,oC,x,fc,cc,0);
				}
			}
		}
	}
}
