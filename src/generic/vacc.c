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
// File:     vcell.c
//
// Purpose:  Class Valist: methods.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbs/vacc.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VACC)
#endif /* if !defined(VINLINE_VACC) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_ctor
//
// Purpose:  Construct the accessibility object
// Notes:    probe_radius is the probe radius (in A) for constructing the 
//             solvent-accessible surface
//           nx, ny, nz are the number of cells (in each direction) to divide
//             the system into for faster access
//           nsph is the number of points (on the surface of the sphere) used
//             to assess solvent accessibility
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vacc* Vacc_ctor(Valist *alist, double probe_radius, int nx,
    int ny, int nz, int nsphere) {


    Vacc *thee = VNULL;

    /* Set up the structure */
    thee = Vram_ctor( 1, sizeof(Vacc) );
    VASSERT( thee != VNULL);
    VASSERT( Vacc_ctor2(thee, alist, probe_radius, nx, ny, nz, nsphere));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_ctor2
//
// Purpose:  Construct the accessibility table object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vacc_ctor2(Vacc *thee, Valist *alist, double probe_radius,
    int nx, int ny, int nz, int nsphere) {

    /* Grid variables */
    int i;
    double x, y, z;
    double rx, ry, rz;
    double x_max, y_max, z_max;
    double x_min, y_min, z_min;
    int ii, jj, kk;
    int i_min, j_min, k_min;
    int i_max, j_max, k_max;
    /* Natural grid coordinate (array position) */
    int ui;
    /* Atom radius */
    double max_radius;
    double tot_r;
    /* A list of how many atoms we've already stored at each grid point */
    int *nputatoms;
    /* Atom list */
    Vatom *atom_list;

    /* Set up grid dimensions */
    thee->nx = nx;
    thee->ny = ny;
    thee->nz = ny;
    thee->nsphere = nsphere;
    thee->n = nx*ny*nz;

    /* Allocate space */
    if ((thee->natoms = Vram_ctor(thee->n,(sizeof(int)))) == VNULL) {
        fprintf(stderr, "Vacc_ctor2: Failed to allocate space.\n");
        return 0;
    }
    for (i=0; i<thee->n; i++) (thee->natoms)[i] = 0;

    if ((thee->atoms = Vram_ctor(thee->n,(sizeof(Vatom **)))) == VNULL) {
        fprintf(stderr, "Vacc_ctor2: Failed to allocate space.\n");
        return 0;
    }
    for (i=0; i<thee->n; i++) (thee->atoms)[i] = VNULL;

    /* Get atom list */
    VASSERT(alist != VNULL);
    atom_list = Valist_getAtomList(alist);

    /* Find dimensions of protein and atoms*/
    x_max = y_max = z_max = -VLARGE;
    x_min = y_min = z_min = VLARGE;
    max_radius = -1.0;
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        x = (Vatom_getPosition(atom_list + i))[0];
        y = (Vatom_getPosition(atom_list + i))[1];
        z = (Vatom_getPosition(atom_list + i))[2];
        if (x < x_min) x_min = x;
        if (y < y_min) y_min = y;
        if (z < z_min) z_min = z;
        if (x > x_max) x_max = x;
        if (y > y_max) y_max = y;
        if (z > z_max) z_max = z;
        if (Vatom_getRadius(atom_list + i) > max_radius) {
            max_radius = Vatom_getRadius(atom_list + i );
        }
    }

    /* Set the probe radius */
    thee->probe_radius = probe_radius;

    /* Set up grid spacings, 2.84 > 2*sqrt(2) */
    thee->hx = (x_max - x_min + 2.84*(max_radius+thee->probe_radius))/(thee->nx-2);
    thee->hy = (y_max - y_min + 2.84*(max_radius+thee->probe_radius))/(thee->ny-2);
    thee->hz = (z_max - z_min + 2.84*(max_radius+thee->probe_radius))/(thee->nz-2);
 
    /* Inflate the grid a bit 1.42 > sqrt(2) */
    (thee->grid_lower_corner)[0] = x_min-1.42*(max_radius+thee->probe_radius)-thee->hx;
    (thee->grid_lower_corner)[1] = y_min-1.42*(max_radius+thee->probe_radius)-thee->hy;
    (thee->grid_lower_corner)[2] = z_min-1.42*(max_radius+thee->probe_radius)-thee->hz;

    /* Find out how many atoms are associated with each grid point */
    for (i=0;i<Valist_getNumberAtoms(alist);i++) { 
 
        /* Get the position in the grid's frame of reference */
        x = (Vatom_getPosition(atom_list+i))[0] 
                - (thee->grid_lower_corner)[0];
        y = (Vatom_getPosition(atom_list+i))[1] 
                - (thee->grid_lower_corner)[1];
        z = (Vatom_getPosition(atom_list+i))[2] 
                - (thee->grid_lower_corner)[2];

        /* Get the range the atom radius + probe radius spans */
        tot_r = Vatom_getRadius(&(atom_list[i])) + thee->probe_radius;
    
        /* Calculate the range of grid points the inflated atom spans in the x 
           direction .
           I admit, this method for determing the accessible grid points isn't
           very fancy, but I kept getting crap for answers with Pythagorean's
           Thm */
        i_max = (int)( ceil( (x + tot_r)/(thee->hx) ));
        i_min = (int)(floor( (x - tot_r)/(thee->hx) ));
        j_max = (int)( ceil( (y + tot_r)/(thee->hy) ));
        j_min = (int)(floor( (y - tot_r)/(thee->hy) ));
        k_max = (int)( ceil( (z + tot_r)/(thee->hz) ));
        k_min = (int)(floor( (z - tot_r)/(thee->hz) ));
 
        /* Now find and assign the grid points */
        for ( ii = i_min; ii <= i_max; ii++) {
            for ( jj = j_min; jj <= j_max; jj++) {
                for ( kk = k_min; kk <= k_max; kk++) {
                    rx = VABS((double)(ii)*(thee->hx) - x);
                    if (rx > thee->hx) rx = rx - thee->hx;
                    ry = VABS((double)(jj)*(thee->hy) - y);
                    if (ry > thee->hy) ry = ry - thee->hy;
                    rz = VABS((double)(kk)*(thee->hz) - z);
                    if (rz > thee->hz) rz = rz - thee->hz;
                    /* Fudge to correct for numerical problems */
                    if ((rx*rx + ry*ry + rz*rz) <= (tot_r*tot_r )) {
                        ui = (thee->nz)*(thee->ny)*ii + (thee->nz)*jj + kk;
                       (thee->natoms[ui])++;
                    }  /* if ... */
                } /* for (kk ... */
            } /* for (jj ... */
        } /* for (ii ... */
    } /* for i =0;i<Valist_getNumberAtoms(alist);i++) */

    /* Allocate the space to store the pointers to the atoms */
    for (i=0;i<thee->n;i++) {
        if ((thee->natoms)[i] > 0) {
            thee->atoms[i] = Vram_ctor((thee->natoms)[i],sizeof(Vatom *));
            VASSERT(thee->atoms[i] != VNULL);
        }
    }

    if (( nputatoms = Vram_ctor(thee->n,sizeof(int))) == VNULL) {
        fprintf(stderr,"Vacc_ctor2: Failed to allocate space.\n");
        return 0;
    }
    for (i=0;i<thee->n;i++) nputatoms[i] = 0;
    for (i=0;i<alist->number;i++) { 
 
        /* Get the position in the grid's frame of reference */
        x = (Vatom_getPosition(atom_list +i))[0] - (thee->grid_lower_corner)[0];
        y = (Vatom_getPosition(atom_list +i))[1] - (thee->grid_lower_corner)[1];
        z = (Vatom_getPosition(atom_list +i))[2] - (thee->grid_lower_corner)[2];

    /* Get the range the atom radius + probe radius spans */
    tot_r = ((alist->atoms)[i]).radius + thee->probe_radius;
    
    /* Calculate the range of grid points the inflated atom spans in the x 
       direction 
       
       I admit, this method for determing the accessible grid points isn't
       very fancy, but I kept getting crap for answers with Pythagorean's
       Thm */

    i_max = (int)( ceil( (x + tot_r)/(thee->hx) ));
    i_min = (int)(floor( (x - tot_r)/(thee->hx) ));
    j_max = (int)( ceil( (y + tot_r)/(thee->hy) ));
    j_min = (int)(floor( (y - tot_r)/(thee->hy) ));
    k_max = (int)( ceil( (z + tot_r)/(thee->hz) ));
    k_min = (int)(floor( (z - tot_r)/(thee->hz) ));

    /* Now find and assign the grid points */
    for ( ii = i_min; ii <= i_max; ii++) {
      for ( jj = j_min; jj <= j_max; jj++) {
        for ( kk = k_min; kk <= k_max; kk++) {
          rx = VABS((double)(ii)*(thee->hx) - x);
          if (rx > thee->hx) rx = rx - thee->hx;
          ry = VABS((double)(jj)*(thee->hy) - y);
          if (ry > thee->hy) ry = ry - thee->hy;
          rz = VABS((double)(kk)*(thee->hz) - z);
          if (rz > thee->hz) rz = rz - thee->hz;
          /* Fudge to correct for numerical problems */
          if ((rx*rx + ry*ry + rz*rz) <= (tot_r*tot_r )) {
            ui = (thee->nz)*(thee->ny)*ii + (thee->nz)*jj + kk;
            thee->atoms[ui][nputatoms[ui]] = &(alist->atoms[i]);
            (nputatoms[ui])++;
          }                                                
        }
      }
    }
  } 

  Vram_dtor((Vram **)&nputatoms,thee->n,sizeof(int));

  /* Set up points on the surface of a sphere */
  thee->sphere = Vacc_sphere(thee, &(thee->nsphere), thee->probe_radius);
  VASSERT(thee->sphere != VNULL);

  return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_dtor
//
// Purpose:  Clean up
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vacc_dtor(Vacc **thee) {
    
    if ((*thee) != VNULL) {
        Vacc_dtor2(*thee);
        Vram_dtor((Vram **)thee, 1, sizeof(Vacc));
        (*thee) = VNULL;
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_dtor2
//
// Purpose:  Clean up
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vacc_dtor2(Vacc *thee) {

    int i;
    for (i=0; i<thee->n; i++) {
        Vram_dtor((Vram **)&((thee->atoms)[i]),(thee->natoms)[i],sizeof(Vatom *));
    }
    Vram_dtor((Vram **)&(thee->natoms),thee->n, sizeof(int));
    Vram_dtor((Vram **)&(thee->sphere),thee->nsphere, sizeof(Vec3));
}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_vdwAcc
//
// Purpose:  Determines if a point is within the union of the atomic spheres
//           (with radii equal to their van der Waals radii).
//           Returns 1 if accessible (outside the molecular volume).
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vacc_vdwAcc(Vacc *thee, Vec3 center) {

    int centeri, centerj, centerk;  /* Grid-based coordinates */
    int ui;                         /* Natural array coordinates */
    int iatom;                      /* Counters */
    double dist;
    Vec3 vec;

    /* Convert to grid based coordinates */
    centeri = VRINT( (center[0] - (thee->grid_lower_corner)[0])/thee->hx);
    centerj = VRINT( (center[1] - (thee->grid_lower_corner)[1])/thee->hy);
    centerk = VRINT( (center[2] - (thee->grid_lower_corner)[2])/thee->hz);
   
    /* Check to make sure we're on the grid; if not, we're definitely 
     * accessible */ 
    if ((centeri < 0) || (centeri >= thee->nx) || 
        (centerj < 0) || (centerj >= thee->ny) || 
        (centerk < 0) || (centerk >= thee->nz)) return 1;
   
    /* If we're still here, then we need to check each atom until we find an
     * overlap at which point we can determine that the point is not 
     * accessible */
    ui = (thee->nz)*(thee->ny)*centeri + (thee->nz)*centerj + centerk;
    for (iatom=0;iatom<(thee->natoms)[ui];iatom++) {
        vec[0] = (Vatom_getPosition((thee->atoms)[ui][iatom]))[0];
        vec[1] = (Vatom_getPosition((thee->atoms)[ui][iatom]))[1];
        vec[2] = (Vatom_getPosition((thee->atoms)[ui][iatom]))[2];
        dist = Vec3_dif2(center,vec);
        if (dist < Vatom_getRadius((thee->atoms)[ui][iatom])) return 0;
    }

    /* If we're still here, then the point is accessible */
    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_ivdwAcc
//
// Purpose:  Determines if a point is within the union of the spheres centered
//           at the atomic centers with radii equal to the sum of their van 
//           der Waals radii and the probe radius.
//           Returns 1 if accessible (outside the molecular volume).
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vacc_ivdwAcc(Vacc *thee, Vec3 center) {

    int centeri, centerj, centerk;  /* Grid-based coordinates */
    int ui;                         /* Natural array coordinates */
    int iatom;                      /* Counters */
    double dist;
    Vec3 vec;

    /* Convert to grid based coordinates */
    centeri = VRINT( (center[0] - (thee->grid_lower_corner)[0])/thee->hx);
    centerj = VRINT( (center[1] - (thee->grid_lower_corner)[1])/thee->hy);
    centerk = VRINT( (center[2] - (thee->grid_lower_corner)[2])/thee->hz);

    /* Check to make sure we're on the grid; if not, we're definitely
     * accessible */
    if ((centeri < 0) || (centeri >= thee->nx) || \
        (centerj < 0) || (centerj >= thee->ny) || \
        (centerk < 0) || (centerk >= thee->nz)) return 1;

    /* If we're still here, then we need to check each atom until we find an
     * overlap at which point we can determine that the point is not
     * accessible */
    ui = (thee->nz)*(thee->ny)*centeri + (thee->nz)*centerj + centerk;
    for (iatom=0;iatom<(thee->natoms)[ui];iatom++) {
        vec[0] = (Vatom_getPosition((thee->atoms)[ui][iatom]))[0];
        vec[1] = (Vatom_getPosition((thee->atoms)[ui][iatom]))[1];
        vec[2] = (Vatom_getPosition((thee->atoms)[ui][iatom]))[2];
        dist = Vec3_dif2(center,vec);
        if (dist < 
           (Vatom_getRadius((thee->atoms)[ui][iatom])+thee->probe_radius) )
          return 0;
    }

    /* If we're still here, then the point is accessible */
    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_molAcc
//
// Purpose:  Determine accessibility of a probe (of radius probe_radius)
//           at a given point, given a collection of atomic spheres.
//           Returns 1 if accessible.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vacc_molAcc(Vacc *thee, Vec3 center) {

    int ipt;
    Vec3 vec;

    /* ******* CHECK IF OUTSIDE ATOM+PROBE RADIUS SURFACE ***** */
    if (Vacc_ivdwAcc(thee, center)) return 1;

    /* ******* CHECK IF INSIDE ATOM RADIUS SURFACE ***** */
    if (!Vacc_vdwAcc(thee, center)) return 0;

    /* ******* CHECK IF OUTSIDE MOLECULAR SURFACE ***** */
    /* Let S be the sphere of radius probe_radius centered at the point we are
     * testing.  We are outside the molecular surface if there is a point on
     * the surface of S that is outside the atom+probe radius surface */
    VASSERT(thee->sphere != VNULL);
    for (ipt=0; ipt<thee->nsphere; ipt++) {
        vec[0] = thee->sphere[ipt][0] + center[0];
        vec[1] = thee->sphere[ipt][1] + center[1];
        vec[2] = thee->sphere[ipt][2] + center[2];
        if (Vacc_ivdwAcc(thee, vec)) return 1;
    }

    /* If all else failed, we are not inside the molecular surface */
    return 0;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_sphere
//
// Purpose:  Generates thee->npts somewhat uniformly distributed across a
//           sphere of specified radius centered at the origin.  Returns a
//           (npts x 3) double array, which you are resposible for destroying,
//           of approximatel the specified number of points; the actual number
//           is stored in the argument npts.  This routine was shamelessly
//           ripped off of sphere.F from UHBD as developed by Michael K. Gilson
//
// Author:   Nathan Baker (original FORTRAN routine from UHBD by Michael
//           Gilson)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double** Vacc_sphere(Vacc *thee, int *npts, double radius) {

    double **points = VNULL;
    int nactual, i, itheta, ntheta, iphi, nphimax, nphi;
    double frac;
    double sintheta, costheta, theta, dtheta;
    double sinphi, cosphi, phi, dphi;

    frac = ((double)(*npts))/4.0;
    ntheta = VRINT(VSQRT(Vunit_pi*frac));
    dtheta = Vunit_pi/((double)(ntheta));
    nphimax = 2*ntheta;

    /* COUNT THE ACTUAL NUMBER OF POINTS TO BE USED */
    nactual = 0;
    for (itheta=0; itheta<ntheta; itheta++) {
        theta = dtheta*((double)(itheta));
        sintheta = VSIN(theta);
        costheta = VCOS(theta);
        nphi = VRINT(sintheta*nphimax);
        nactual += nphi;
    }

    /* ALLOCATE THE SPACE FOR THE POINTS */
    points = Vram_ctor(nactual, sizeof(double *));
    VASSERT(points != VNULL);
    for (i=0; i<nactual; i++) {
        points[i] = Vram_ctor(3, sizeof(double));
        VASSERT(points[i] != VNULL);
    }

    /* ASSIGN THE POINTS */
    nactual = 0;
    for (itheta=0; itheta<ntheta; itheta++) {
        theta = dtheta*((double)(itheta));
        sintheta = VSIN(theta);
        costheta = VCOS(theta);
        nphi = VRINT(sintheta*nphimax);
        if (nphi != 0) {
            dphi = 2*Vunit_pi/((double)(nphi));
            for (iphi=0; iphi<nphi; iphi++) {
                phi = dphi*((double)(iphi));
                sinphi = VSIN(phi);
                cosphi = VCOS(phi);
                points[nactual][0] = radius * cosphi * sintheta;
                points[nactual][1] = radius * sinphi * sintheta;
                points[nactual][2] = radius * costheta;
                nactual++;
            }
        }
    }

    *npts = nactual;
    return points;
}
