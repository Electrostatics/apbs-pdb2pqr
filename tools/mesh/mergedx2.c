/**
*  @file     mergedx2.c
 *  @author  David Gohara, Stephen Bond and Nathan Baker
 *  @brief   Program that merges OpenDX files
 *  @version $Id: mergedx.c 1033 2007-02-25 17:08:22Z sdg0919 $
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
 * Copyright (c) 2002-2007.  Washington University in St. Louis.
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

#include <unistd.h>

#include "apbscfg.h"
#include "apbs/apbs.h"

#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define INTERVAL(x,a,b) (((x) >= (a)) && ((x) <= (b)))

VEMBED(rcsid="$Id: mergedx.c 1033 2006-12-29 17:08:22Z sdg0919 $")

VPRIVATE int Vgrid_readDXhead(Vgrid *thee,
  const char *iodev, const char *iofmt, const char *thost, const char *fname);
VPRIVATE int Vgrid_value2(Vgrid *thee, double pt[3], double *value);

VPRIVATE char *MCwhiteChars = " =,;\t\n";
VPRIVATE char *MCcommChars  = "#%";

void usage(){
	
	Vnm_print(1,"mergedx2 [FLAGS] file1.dx [file2.dx ...]\n"
				"-o		Output file					(default: gridmerged.dx)\n"
				"-r		Resolution of gridpoints	(default: 1.0 Angstroms)\n"
				"-b		Bounds of output map as: xmin ymin zmin xmax ymax zmax\n"
				"									(default: calculates full map)\n"
				"-s		Print bounds of merged input dx files. Doesn't generate a merged map.\n"
				"									(-s is exclusive of the other FLAGS)\n"
				"-h		Print this message\n"
		   );
	
}

int main(int argc, char **argv) {

    /* *************** VARIABLES ******************* */
    int i, j, k,spec,warn;
    int nx, ny, nz, count, numfnams;
	
    double pt[3],value, res;
	
    double xmin, ymin, zmin;
	double xmax, ymax, zmax;
	
	double xminb, yminb, zminb;
	double xmaxb, ymaxb, zmaxb;
	
	/* We will cache the file names by address. So it needs to be 64-bit clean */
    intptr_t fnams[1024];
	short *carray = VNULL;
	
    char *snam = "# main:  ";
    char outname[80];
	
    Vgrid *grid, *mgrid;
	
	/* Set the default values */
	spec = 0;
	warn = 0;
	res = 1.0;
	xmin = ymin = zmin = 0.0;
	xmax = ymax = zmax = 0.0;
	Vnm_print(1,outname,"gridmerged.dx");
	
	/* Begin processing command line options */
	int ch, ind;
	extern char *optarg;
	extern int optind, opterr, optopt;
	
	/* Check the invocation */
	if(argc <= 1){ usage(); return 1; }
	
	while ((ch = getopt(argc, argv, "r:b:o:s:h")) != -1) {
		switch (ch) {
			case 'r':
				res = atof(optarg);
				break;
			case 'b':
				ind = optind - 1;
				xmin = atof(argv[ind++]);
				ymin = atof(argv[ind++]);
				zmin = atof(argv[ind++]);
				
				xmax = atof(argv[ind++]);
				ymax = atof(argv[ind++]);
				zmax = atof(argv[ind++]);
				
				optind = ind;
				break;
			case 'o':
				strcpy(outname,optarg);
				break;
			case 's':
				spec = 1;
				break;
			case 'h':
				usage();
				return 0;
				break;
			default:
				break;
		}
	}
	
	numfnams = 0;
	if (optind < argc) {
		while (optind < argc){
			fnams[numfnams] = *(&argv[optind++]);
			numfnams += 1;
		}
	}
	
	/* Start the I/O processing */
	Vio_start();
	
	/* For now we only allow one resolution on all three axes */
	double resx = res;
	double resy = res;
	double resz = res;
	
    /* *********** PREPARE MERGED GRID OBJECT ************* */
    mgrid = Vgrid_ctor(nx, ny, nz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
    mgrid->xmin = VLARGE; mgrid->xmax = -VLARGE;
    mgrid->ymin = VLARGE; mgrid->ymax = -VLARGE;
    mgrid->zmin = VLARGE; mgrid->zmax = -VLARGE;
	
    /* *************** GET FILE HEADERS ******************* */
    Vnm_print(1, "%s Reading Headers...\n",snam);
    grid = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
    for(count=0; count<numfnams; count++) {
        Vnm_print(0, "%s  Reading header from %s...\n",snam, fnams[count]);
        Vgrid_readDXhead(grid, "FILE", "ASC", VNULL, fnams[count]);
		
        /* set the merged grid bounds to include all the subgrids */
        if( grid->xmin < mgrid->xmin ) mgrid->xmin = grid->xmin;
        if( grid->xmax > mgrid->xmax ) mgrid->xmax = grid->xmax;
        if( grid->ymin < mgrid->ymin ) mgrid->ymin = grid->ymin;
        if( grid->ymax > mgrid->ymax ) mgrid->ymax = grid->ymax;
        if( grid->zmin < mgrid->zmin ) mgrid->zmin = grid->zmin;
        if( grid->zmax > mgrid->zmax ) mgrid->zmax = grid->zmax;
		
		if( grid->hx > res || grid->hy > res || grid->hzed > res ) warn = 1;
    }
	
	if(warn){
		Vnm_print(1,"WARNING: The specified output resolution is greater than the\n"
					"		 resolution of the input DX files. Upsampling.......\n");
	}
	
	/* Cache the bounds for comparison later */
	xminb = mgrid->xmin; yminb = mgrid->ymin; zminb = mgrid->zmin;
	xmaxb = mgrid->xmax; ymaxb = mgrid->ymax; zmaxb = mgrid->zmax;
	
	/* Adjust the boundaries of the grid to any specified by the user */
	if(xmin != 0.0) mgrid->xmin = xmin;
	if(ymin != 0.0) mgrid->ymin = ymin;
	if(zmin != 0.0) mgrid->zmin = zmin;
	
	if(xmax != 0.0) mgrid->xmax = xmax;
	if(ymax != 0.0) mgrid->ymax = ymax;
	if(zmax != 0.0) mgrid->zmax = zmax;
	
	/* Now check the boundaries the user specified (if any) 
		to make sure they fit within the original
	 */
	if((mgrid->xmin < xminb) ||
	   (mgrid->ymin < yminb) ||
	   (mgrid->zmin < zminb) ||
	   (mgrid->xmax > xmaxb) ||
	   (mgrid->ymax > ymaxb) ||
	   (mgrid->zmax > zmaxb))
	{
		Vnm_print(1,"\nError: The bounds requested do not fall within the bounds of the specified grid\n"
					"You specified <xmin> <ymin> <zmin> <xmax> <ymax> <zmax>: %lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
					"The input DX files provided                            : %lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
					mgrid->xmin,mgrid->ymin,mgrid->zmin,mgrid->xmax,mgrid->ymax,mgrid->zmax,
					xminb,yminb,zminb,xmaxb,ymaxb,zmaxb
				  );
		return 1;
	}
	   
    /* set the grid increment for the merged grid */
	mgrid->nx	= VFLOOR(((mgrid->xmax - mgrid->xmin) / resx) + 1.5);
	mgrid->ny	= VFLOOR(((mgrid->ymax - mgrid->ymin) / resy) + 1.5);
	mgrid->nz	= VFLOOR(((mgrid->zmax - mgrid->zmin) / resz) + 1.5);
	
	mgrid->hx   = (mgrid->xmax - mgrid->xmin) / (mgrid->nx-1);
	mgrid->hy   = (mgrid->ymax - mgrid->ymin) / (mgrid->ny-1);
	mgrid->hzed = (mgrid->zmax - mgrid->zmin) / (mgrid->nz-1);
	
    /* print out the dimensions of the merged grid */
    Vnm_print(1, "%s Dimensions of the merged grid\n",snam);
    Vnm_print(1, "%s nx = %d, ny = %d, nz = %d\n",snam,
              mgrid->nx, mgrid->ny, mgrid->nz);
    Vnm_print(1, "%s hx = %lf, hy = %lf, hz = %lf\n",snam,
              mgrid->hx, mgrid->hy, mgrid->hzed);
    Vnm_print(1, "%s xmin = %lf, ymin = %lf, zmin = %lf\n",snam,
              mgrid->xmin, mgrid-> ymin, mgrid->zmin);
    Vnm_print(1, "%s xmax = %lf, ymax = %lf, zmax = %lf\n",snam,
              mgrid->xmax, mgrid-> ymax, mgrid->zmax);
	
	if(spec) return 0;
	
    mgrid->data = (double *)
       Vmem_malloc(mgrid->mem,(mgrid->nx*mgrid->ny*mgrid->nz),sizeof(double));
    mgrid->ctordata = 1;
	
    carray = (short *)
       Vmem_malloc(VNULL, (mgrid->nx*mgrid->ny*mgrid->nz), sizeof(short) );
	
    /* initialize the data of the merged grid with zeros */
    nx = mgrid->nx;
    ny = mgrid->ny;
    nz = mgrid->nz;
    for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++) {
            for (k=0; k<nz; k++) {
                (mgrid->data)[IJK(i,j,k)] = 0.0;
                carray[IJK(i,j,k)] = 0;
            }
        }
    }

    /* ************** MERGE THE GRID FILES **************** */
    Vnm_print(1, "%s Reading and Merging...\n",snam);
    for (count=0; count<numfnams; count++) {
        Vgrid_readDX(grid, "FILE", "ASC", VNULL, fnams[count]);
        xmin = grid->xmin - grid->hx   - VSMALL;
        ymin = grid->ymin - grid->hy   - VSMALL;
        zmin = grid->zmin - grid->hzed - VSMALL;
        xmax = grid->xmax + grid->hx   + VSMALL;
        ymax = grid->ymax + grid->hy   + VSMALL;
        zmax = grid->zmax + grid->hzed + VSMALL;
        for (i=0; i<nx; i++) {
            pt[0] = mgrid->xmin + i*mgrid->hx;
            if(INTERVAL(pt[0],xmin,xmax)) {
                for (j=0; j<ny; j++) {
                    pt[1] = mgrid->ymin + j*mgrid->hy;
                    if(INTERVAL(pt[1],ymin,ymax)) {
                        for (k=0; k<nz; k++) {
                            pt[2] = mgrid->zmin + k*mgrid->hzed;
                            if(INTERVAL(pt[2],zmin,zmax)) {
                                if (Vgrid_value2(grid, pt, &value)) {
                                    (mgrid->data)[IJK(i,j,k)] += value;
                                    carray[IJK(i,j,k)] += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        Vmem_free(grid->mem,(grid->nx*grid->ny*grid->nz), sizeof(double),
                  (void **)&(grid->data));
        grid->readdata = 0;
        grid->ctordata = 0;
    }
    Vgrid_dtor( &grid );

    mgrid->readdata = 1;

    /* Check for skipped points, and account for overlap */
    nx = mgrid->nx;
    ny = mgrid->ny;
    nz = mgrid->nz;
    for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++) {
            for (k=0; k<nz; k++) {
                if ( carray[IJK(i,j,k)] >= 1 ) {
                    (mgrid->data)[IJK(i,j,k)] /= carray[IJK(i,j,k)];
                } else {
                    Vnm_print(2,"%s Warning: Gap in subgrids at point (%g,%g,%g)\n",
							  snam,
                              mgrid->xmin + i*mgrid->hx,
                              mgrid->ymin + j*mgrid->hy,
                              mgrid->zmin + k*mgrid->hzed );
                }
            }
        }
    }

    /* ************** WRITE THE MERGED GRID **************** */
    Vnm_print(1, "%s Writing...\n",snam);
    Vnm_print(0, "%s  Writing merged data to %s...\n",snam,outname);
    Vgrid_writeDX(mgrid, "FILE", "ASC", VNULL, outname,"mergedx",VNULL);

    Vmem_free(VNULL,(mgrid->nx*mgrid->ny*mgrid->nz), sizeof(short),
              (void **)&(carray));
    Vgrid_dtor( &mgrid );

    return 0;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vgrid_value2
//
// Author:   Nathan Baker and Stephen Bond
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vgrid_value2(Vgrid *thee, double pt[3], double *value) {

    int nx, ny, nz, ihi, jhi, khi, ilo, jlo, klo;
    double hx, hy, hzed, xmin, ymin, zmin, ifloat, jfloat, kfloat;
    double u, dx, dy, dz;

    VASSERT(thee != VNULL);
    VASSERT(thee->ctordata || thee->readdata);

    nx = thee->nx;
    ny = thee->ny;
    nz = thee->nz;
    hx = thee->hx;
    hy = thee->hy;
    hzed = thee->hzed;
    xmin = thee->xmin;
    ymin = thee->ymin;
    zmin = thee->zmin;

    u = 0;

    ifloat = (pt[0] - xmin)/hx;
    jfloat = (pt[1] - ymin)/hy;
    kfloat = (pt[2] - zmin)/hzed;
    ihi = (int)ceil(ifloat);
    jhi = (int)ceil(jfloat);
    khi = (int)ceil(kfloat);
    ilo = (int)floor(ifloat);
    jlo = (int)floor(jfloat);
    klo = (int)floor(kfloat);

    /* If the point is outside the mesh, push it to the mesh */
    if ( ilo < 0 ) {
        ilo = 0;
        ihi = ilo + 1;
        ifloat = (double)(ilo);
    } else if ( ihi >= nx ) {
        ihi = nx - 1;
        ilo = ihi - 1;
        ifloat = (double)(ihi);
    }
    if ( jlo < 0 ) {
        jlo = 0;
        jhi = jlo + 1;
        jfloat = (double)(jlo);
    } else if ( jhi >= ny ) {
        jhi = ny - 1;
        jlo = jhi - 1;
        jfloat = (double)(jhi);
    }
    if ( klo < 0 ) {
        klo = 0;
        khi = klo + 1;
        kfloat = (double)(klo);
    } else if ( khi >= nz ) {
        khi = nz - 1;
        klo = khi - 1;
        kfloat = (double)(khi);
    } 

    /* See if we're on the mesh */
    if ((ihi<nx) && (jhi<ny) && (khi<nz) &&
        (ilo>=0) && (jlo>=0) && (klo>=0)) {
        dx = ifloat - (double)(ilo);
        dy = jfloat - (double)(jlo);
        dz = kfloat - (double)(klo);
        u = dx      *dy      *dz      *(thee->data[IJK(ihi,jhi,khi)])
          + dx      *(1.0-dy)*dz      *(thee->data[IJK(ihi,jlo,khi)])
          + dx      *dy      *(1.0-dz)*(thee->data[IJK(ihi,jhi,klo)])
          + dx      *(1.0-dy)*(1.0-dz)*(thee->data[IJK(ihi,jlo,klo)])
          + (1.0-dx)*dy      *dz      *(thee->data[IJK(ilo,jhi,khi)])
          + (1.0-dx)*(1.0-dy)*dz      *(thee->data[IJK(ilo,jlo,khi)])
          + (1.0-dx)*dy      *(1.0-dz)*(thee->data[IJK(ilo,jhi,klo)])
          + (1.0-dx)*(1.0-dy)*(1.0-dz)*(thee->data[IJK(ilo,jlo,klo)]);

        *value = u;
        return 1;

    }

    *value = 0.0;
    return 0;
}

/* ///////////////////////////////////////////////////////////////////////////
   // Routine:  Vgrid_readDXhead
   //
   // Author:   Nathan Baker and Stephen Bond
   /////////////////////////////////////////////////////////////////////////// */
VPRIVATE int Vgrid_readDXhead(Vgrid *thee,
							  const char *iodev, const char *iofmt, const char *thost, const char *fname) {
	
    int itmp;
    double dtmp;
    char tok[VMAX_BUFSIZE];
    char *snam = "Vgrid_readDXhead:";
    Vio *sock;
	
    /* Check to see if the existing data is null and, if not, clear it out */
    if (thee->data != VNULL) {
        Vnm_print(1, "%s  destroying existing data!\n",snam);
        Vmem_free(thee->mem, (thee->nx*thee->ny*thee->nz), sizeof(double),
				  (void **)&(thee->data)); }
    thee->readdata = 0;
    thee->ctordata = 0;
	
    /* Set up the virtual socket */
    sock = Vio_ctor(iodev,iofmt,thost,fname,"r");
    if (sock == VNULL) {
        Vnm_print(2, "%s Problem opening virtual socket %s\n",snam,fname);
        return 0;
    }
    if (Vio_accept(sock, 0) < 0) {
        Vnm_print(2, "%s Problem accepting virtual socket %s\n",snam,fname);
        return 0;
    }
	
    Vio_setWhiteChars(sock, MCwhiteChars);
    Vio_setCommChars(sock, MCcommChars);
	
    /* Read in the DX regular positions */
    /* Get "object" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "object"));
    /* Get "1" */
    VJMPERR2(1 == Vio_scanf(sock, "%d", &itmp));
    /* Get "class" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "class"));
    /* Get "gridpositions" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "gridpositions"));
    /* Get "counts" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "counts"));
    /* Get nx */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%d", &(thee->nx)));
    /* Get ny */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%d", &(thee->ny)));
    /* Get nz */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%d", &(thee->nz)));
    Vnm_print(0, "%s  Grid dimensions %d x %d x %d grid\n",snam,
              thee->nx, thee->ny, thee->nz);
    /* Get "origin" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "origin"));
    /* Get xmin */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &(thee->xmin)));
    /* Get ymin */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &(thee->ymin)));
    /* Get zmin */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &(thee->zmin)));
    Vnm_print(0, "%s  Grid origin = (%g, %g, %g)\n",snam,
              thee->xmin, thee->ymin, thee->zmin);
    /* Get "delta" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "delta"));
    /* Get hx */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &(thee->hx)));
    /* Get 0.0 */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
    VJMPERR1(dtmp == 0.0);
    /* Get 0.0 */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
    VJMPERR1(dtmp == 0.0);
    /* Get "delta" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "delta"));
    /* Get 0.0 */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
    VJMPERR1(dtmp == 0.0);
    /* Get hy */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &(thee->hy)));
    /* Get 0.0 */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
    VJMPERR1(dtmp == 0.0);
    /* Get "delta" */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(!strcmp(tok, "delta"));
    /* Get 0.0 */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
    VJMPERR1(dtmp == 0.0);
    /* Get 0.0 */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &dtmp));
    VJMPERR1(dtmp == 0.0);
    /* Get hz */
    VJMPERR2(1 == Vio_scanf(sock, "%s", tok));
    VJMPERR1(1 == sscanf(tok, "%lf", &(thee->hzed)));
    Vnm_print(0, "%s  Grid spacings = (%g, %g, %g)\n",snam,
              thee->hx, thee->hy, thee->hzed);
    /* calculate grid maxima */
    thee->xmax = thee->xmin + (thee->nx-1)*thee->hx;
    thee->ymax = thee->ymin + (thee->ny-1)*thee->hy;
    thee->zmax = thee->zmin + (thee->nz-1)*thee->hzed;
	
    /* Close off the socket */
    Vio_acceptFree(sock);
    Vio_dtor(&sock);
	
    return 1;
	
VERROR1:
	Vio_dtor(&sock);
    Vnm_print(2, "%s  Format problem with input file <%s>\n",snam,fname);
    return 0;
	
VERROR2:
		Vio_dtor(&sock);
    Vnm_print(2, "%s  I/O problem with input file <%s>\n",snam,fname);
    return 0;
}
