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
// File:     vacc.c
//
// Purpose:  Class Vacc: methods.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "apbs/vacc.h"

#if defined(HAVE_FETK_H)
#include "mc/mc.h"
#endif

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VACC)

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_memChk
//
// Purpose:  Returns the number of bytes used by the specified object.
//
// Author:   Nathan Baker 
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vacc_memChk(Vacc *thee) {
    if (thee == VNULL) return 0;
    return Vmem_bytes(thee->vmem);
}

#endif /* if !defined(VINLINE_VACC) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int ivdwAccExclus(Vacc *thee, double center[3], double radius, int atomID);

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
VPUBLIC Vacc* Vacc_ctor(Valist *alist, double max_radius, int nx,
    int ny, int nz, int nsphere) {


    Vacc *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vacc) );
    VASSERT( thee != VNULL);
    VASSERT( Vacc_ctor2(thee, alist, max_radius, nx, ny, nz, nsphere));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_ctor2
//
// Purpose:  Construct the accessibility table object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vacc_ctor2(Vacc *thee, Valist *alist, double max_radius,
    int nx, int ny, int nz, int nsphere) {

    /* Grid variables */
    int i;
    double x, y, z, *coord, dx, dy, dz, dx2, dy2, dz2;
    double x_max, y_max, z_max;
    double x_min, y_min, z_min;
    int ii, jj, kk, totatoms, dumi;
    int i_min, j_min, k_min;
    int i_max, j_max, k_max;
    /* Natural grid coordinate (array position) */
    int ui;
    /* Atom radius */
    double rmax;
    double rtot, rtot2;
    Vatom *atom;

    VASSERT(alist != VNULL);
    thee->alist = alist;

    /* Set up memory management object */
    thee->vmem = Vmem_ctor("APBS::VACC");

    /* Set up grid dimensions */
    thee->nx = nx;
    thee->ny = ny;
    thee->nz = nz;
    thee->n = nx*ny*nz;
    if ((nx < 3) || (ny < 3) || (nz < 3)) {
        Vnm_print(2, "Vacc_ctor2:  nx, ny, nz must be greater than 2!\n");
        return 0;
    }
    Vnm_print(0, "Vacc_ctor2:  Using %d x %d x %d hash table\n", nx, ny, nz);
 
    /* Set up probe information */
    thee->nsphere = nsphere;
    thee->max_radius = max_radius;
    Vnm_print(0, "Vacc_ctor2:  Constructing sphere...\n");
    thee->sphere = Vacc_sphere(thee, &(thee->nsphere));
    VASSERT(thee->sphere != VNULL);

    /* Allocate space */
    thee->natoms = Vmem_malloc(thee->vmem, thee->n, sizeof(int));
    VASSERT(thee->natoms != VNULL);
    for (i=0; i<thee->n; i++) (thee->natoms)[i] = 0;
    thee->atomIDs = Vmem_malloc(thee->vmem, thee->n, sizeof(int *));
    VASSERT(thee->atomIDs != VNULL);
    for (i=0; i<thee->n; i++) (thee->atomIDs)[i] = VNULL;
    thee->area = Vmem_malloc(thee->vmem, Valist_getNumberAtoms(alist), 
        sizeof(double));
    VASSERT(thee->area != VNULL);
    for (i=0; i<Valist_getNumberAtoms(alist); i++) thee->area[i] = 0;

    /* Find dimensions of protein and atoms */
    x_max = y_max = z_max = -VLARGE;
    x_min = y_min = z_min = VLARGE;
    rmax = -1.0;
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom = Valist_getAtom(alist, i);
        x = (Vatom_getPosition(atom))[0];
        y = (Vatom_getPosition(atom))[1];
        z = (Vatom_getPosition(atom))[2];
        if (x < x_min) x_min = x;
        if (y < y_min) y_min = y;
        if (z < z_min) z_min = z;
        if (x > x_max) x_max = x;
        if (y > y_max) y_max = y;
        if (z > z_max) z_max = z;
        if (Vatom_getRadius(atom) > rmax) rmax = Vatom_getRadius(atom);
    }

    /* Set up grid spacings, 2.84 > 2*sqrt(2) */
    thee->hx = (x_max - x_min + 2.84*(rmax+thee->max_radius))/thee->nx;
    thee->hy = (y_max - y_min + 2.84*(rmax+thee->max_radius))/thee->ny;
    thee->hzed = (z_max - z_min + 2.84*(rmax+thee->max_radius))/thee->nz;
 
    /* Inflate the grid a bit 1.42 > sqrt(2) */
#if 0
    (thee->grid_lower_corner)[0] = x_min-1.42*(rmax+thee->max_radius)-thee->hx;
    (thee->grid_lower_corner)[1] = y_min-1.42*(rmax+thee->max_radius)-thee->hy;
    (thee->grid_lower_corner)[2] = z_min-1.42*(rmax+thee->max_radius)-thee->hzed;
#else
    (thee->grid_lower_corner)[0] = x_min-1.42*(rmax+thee->max_radius);
    (thee->grid_lower_corner)[1] = y_min-1.42*(rmax+thee->max_radius);
    (thee->grid_lower_corner)[2] = z_min-1.42*(rmax+thee->max_radius);
#endif

    /* Find out how many atoms are associated with each grid point */
    totatoms = 0;
    for (i=0;i<Valist_getNumberAtoms(alist);i++) { 
        atom = Valist_getAtom(alist, i);
        /* Get the position in the grid's frame of reference */
        coord = Vatom_getPosition(atom);
        x = coord[0] - (thee->grid_lower_corner)[0];
        y = coord[1] - (thee->grid_lower_corner)[1];
        z = coord[2] - (thee->grid_lower_corner)[2];

        /* Get the range the atom radius + probe radius spans */
        rtot = Vatom_getRadius(atom) + thee->max_radius;
    
        /* Calculate the range of grid points the inflated atom spans in the x 
         * direction. */
        i_max = (int)( ceil( (x + rtot)/(thee->hx) ));
        i_max = VMIN2(i_max, nx-1);
        i_min = (int)(floor( (x - rtot)/(thee->hx) ));
        i_min = VMAX2(i_min, 0);
        j_max = (int)( ceil( (y + rtot)/(thee->hy) ));
        j_max = VMIN2(j_max, ny-1);
        j_min = (int)(floor( (y - rtot)/(thee->hy) ));
        j_min = VMAX2(j_min, 0);
        k_max = (int)( ceil( (z + rtot)/(thee->hzed) ));
        k_max = VMIN2(k_max, nz-1);
        k_min = (int)(floor( (z - rtot)/(thee->hzed) ));
        k_min = VMAX2(k_min, 0);

        Vnm_print(0, "VACC DEBUG: %d <= i <= %d\n", i_min, i_max);
        Vnm_print(0, "VACC DEBUG: %d <= j <= %d\n", j_min, j_max);
        Vnm_print(0, "VACC DEBUG: %d <= k <= %d\n", k_min, k_max);

        /* Now find and assign the grid points */
        for ( ii = i_min; ii <= i_max; ii++) {
            for ( jj = j_min; jj <= j_max; jj++) {
                for ( kk = k_min; kk <= k_max; kk++) {
                    ui = (thee->nz)*(thee->ny)*ii + (thee->nz)*jj + kk;
                    (thee->natoms[ui])++;
                    totatoms += thee->natoms[i];
                } 
            } 
        } 
    } /* for i =0:natoms */
    Vnm_print(0, "Vacc_ctor2:  Have %d atom entries\n", totatoms);

    /* Allocate the space to store the pointers to the atoms */
    for (i=0; i<thee->n; i++) {
        if ((thee->natoms)[i] > 0) {
            thee->atomIDs[i] = Vmem_malloc(thee->vmem, thee->natoms[i],
              sizeof(int));
            VASSERT(thee->atomIDs[i] != VNULL);
        }
        /* Clear the counter for later use */
        thee->natoms[i] = 0;
    }
 
    /* Assign the atoms to grid points */
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom = Valist_getAtom(alist, i);
        /* Get the position in the grid's frame of reference */
        x = (Vatom_getPosition(atom))[0] - (thee->grid_lower_corner)[0];
        y = (Vatom_getPosition(atom))[1] - (thee->grid_lower_corner)[1];
        z = (Vatom_getPosition(atom))[2] - (thee->grid_lower_corner)[2];

        /* Get the range the atom radius + probe radius spans */
        rtot = Vatom_getRadius(atom) + thee->max_radius;
        rtot2 = VSQR(rtot);

        /* Now find and assign the grid points */
        i_max = (int)( ceil( (x + rtot)/(thee->hx) ));
        i_max = VMIN2(i_max, nx-1);
        i_min = (int)(floor( (x - rtot)/(thee->hx) ));
        i_min = VMAX2(i_min, 0);
        j_max = (int)( ceil( (y + rtot)/(thee->hy) ));
        j_max = VMIN2(j_max, ny-1);
        j_min = (int)(floor( (y - rtot)/(thee->hy) ));
        j_min = VMAX2(j_min, 0);
        k_max = (int)( ceil( (z + rtot)/(thee->hzed) ));
        k_max = VMIN2(k_max, nz-1);
        k_min = (int)(floor( (z - rtot)/(thee->hzed) ));
        k_min = VMAX2(k_min, 0);
        /* Now find and assign the grid points */
        for ( ii = i_min; ii <= i_max; ii++) {
            for ( jj = j_min; jj <= j_max; jj++) {
                for ( kk = k_min; kk <= k_max; kk++) {

                    ui = (thee->nz)*(thee->ny)*ii + (thee->nz)*jj + kk;
                    thee->atomIDs[ui][thee->natoms[ui]] = i;
                    (thee->natoms[ui])++;
                }
            }
        }
    } /* for i =0:natoms */


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
        Vmem_free(VNULL, 1, sizeof(Vacc), (void **)thee);
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
	if (thee->natoms[i] > 0)  Vmem_free(thee->vmem, (thee->natoms)[i],
          sizeof(int), (void **)&(thee->atomIDs[i]));
    }
    Vmem_free(thee->vmem, thee->n, sizeof(int *), (void **)&(thee->atomIDs));
    Vmem_free(thee->vmem, thee->n, sizeof(int), (void **)&(thee->natoms));
    for (i=0; i<thee->nsphere; i++) 
      Vmem_free(thee->vmem, 3, sizeof(double), (void **)&(thee->sphere[i]));
    Vmem_free(thee->vmem, thee->nsphere, sizeof(double *), 
      (void **)&(thee->sphere));
    Vmem_free(thee->vmem, Valist_getNumberAtoms(thee->alist),
      sizeof(double), (void **)&(thee->area));

    Vmem_dtor(&(thee->vmem));
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
VPUBLIC int Vacc_vdwAcc(Vacc *thee, double center[3]) {

    int centeri, centerj, centerk;  /* Grid-based coordinates */
    int ui;                         /* Natural array coordinates */
    int iatom;                      /* Counters */
    double dist;
    Vatom *atom;
    double *apos;

    /* Convert to grid based coordinates */
    centeri = (int)( (center[0] - (thee->grid_lower_corner)[0])/thee->hx);
    centerj = (int)( (center[1] - (thee->grid_lower_corner)[1])/thee->hy);
    centerk = (int)( (center[2] - (thee->grid_lower_corner)[2])/thee->hzed);
   
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
        atom = Valist_getAtom(thee->alist, thee->atomIDs[ui][iatom]);
        apos = Vatom_getPosition(atom);
        dist = VSQR(center[0]-apos[0]) + VSQR(center[1]-apos[1])
               + VSQR(center[2]-apos[2]);
        if (dist < VSQR(Vatom_getRadius(atom))) return 0;
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
VPUBLIC int Vacc_ivdwAcc(Vacc *thee, double center[3], double radius) {

    return ivdwAccExclus(thee, center, radius, -1);

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  ivdwAccExclus
//
// Purpose:  Determines if a point is within the union of the spheres centered
//           at the atomic centers with radii equal to the sum of their van
//           der Waals radii and the probe radius.  Does not include
//           contributions from the specified atom.
//
// Args:     center => point to be tested
//           radius => radius to inflate by
//           atomID  => atom to ignore (-1 if none)
//           
//           Returns 1 if accessible (outside the molecular volume).
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int ivdwAccExclus(Vacc *thee, double center[3], double radius, 
  int atomID) {

    int centeri, centerj, centerk;  /* Grid-based coordinates */
    int ui;                         /* Natural array coordinates */
    int iatom;                      /* Counters */
    double dist, *apos;
    Vatom *atom;

    /* We can only test probes with radii less than the max specified */
    VASSERT(thee != VNULL);
    if (radius > thee->max_radius) {
        Vnm_print(2, "Vacc_ivdwAcc: got radius (%g) bigger than max radius (%g)\n", 
          radius, thee->max_radius);
         VASSERT(0);
    }

    /* Convert to grid based coordinates */
    centeri = (int)( (center[0] - (thee->grid_lower_corner)[0])/thee->hx);
    centerj = (int)( (center[1] - (thee->grid_lower_corner)[1])/thee->hy);
    centerk = (int)( (center[2] - (thee->grid_lower_corner)[2])/thee->hzed);

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
        if (thee->atomIDs[ui][iatom] != atomID) {
            atom = Valist_getAtom(thee->alist, thee->atomIDs[ui][iatom]);
            apos = Vatom_getPosition(atom);
            dist = VSQR(apos[0]-center[0]) + VSQR(apos[1]-center[1])
                   + VSQR(apos[2]-center[2]);
            if (dist < VSQR(Vatom_getRadius(atom)+radius)) return 0;
        }
    }

    /* If we're still here, then the point is accessible */
    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_molAcc
//
// Purpose:  Determine accessibility of a probe (of radius radius)
//           at a given point, given a collection of atomic spheres.
//           Returns 1 if accessible.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vacc_molAcc(Vacc *thee, double center[3], double radius) {

    int ipt;
    double vec[3];

    /* ******* CHECK IF OUTSIDE ATOM+PROBE RADIUS SURFACE ***** */
    if (Vacc_ivdwAcc(thee, center, radius)) return 1;

    /* ******* CHECK IF INSIDE ATOM RADIUS SURFACE ***** */
    if (!Vacc_vdwAcc(thee, center)) return 0;

    /* ******* CHECK IF OUTSIDE MOLECULAR SURFACE ***** */
    /* Let S be the sphere of radius radius centered at the point we are
     * testing.  We are outside the molecular surface if there is a point on
     * the surface of S that is outside the atom+probe radius surface */
    VASSERT(thee->sphere != VNULL);
    for (ipt=0; ipt<thee->nsphere; ipt++) {
        vec[0] = radius*thee->sphere[ipt][0] + center[0];
        vec[1] = radius*thee->sphere[ipt][1] + center[1];
        vec[2] = radius*thee->sphere[ipt][2] + center[2];
        if (Vacc_ivdwAcc(thee, vec, radius)) return 1;
    }

    /* If all else failed, we are not inside the molecular surface */
    return 0;
}

#if defined(HAVE_FETK_H)
/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_writeGMV
//
// Purpose:  Write out the chosen accessibility data at each vertex.  The the
//           appropriate isosurface routine, this would generate a
//           representation of the molecular surface as ``seen" by the PBE
//           solver.
//
// Arguments: radius   Radius of sphere to test
//            meth     Plot accessibility for molecular surface (meth=0),
//                     inflated van der Waals (meth=1), or van der Waals
//                     (meth=2)
//            gm       Gem object with mesh data
//            iodev    Device (usually "FILE")
//            iofmt    Format (usually "ASC")
//            iohost   Host   (usually "localhost")
//            iofile   Filename
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vacc_writeGMV(Vacc *thee, double radius, int meth, Gem *gm, 
  char *iodev, char *iofmt, char *iohost, char *iofile) {

    double *accVals[MAXV], coord[3];
    int ivert, icoord;

    for (ivert=0; ivert<MAXV; ivert++) accVals[ivert] = VNULL;
    accVals[0] = (void *)Vmem_malloc(thee->vmem, Gem_numVV(gm), sizeof(double));
    accVals[1] = (void *)Vmem_malloc(thee->vmem, Gem_numVV(gm), sizeof(double));
    for (ivert=0; ivert<Gem_numVV(gm); ivert++) {
        for (icoord=0;icoord<3;icoord++) 
          coord[icoord] = VV_coord(Gem_VV(gm, ivert), icoord);
        if (meth == 0) {
            accVals[0][ivert] = (double)Vacc_molAcc(thee, coord, radius);
            accVals[1][ivert] = (double)Vacc_molAcc(thee, coord, radius);
        } else if (meth == 1) {
            accVals[0][ivert] = (double)Vacc_ivdwAcc(thee, coord, radius);
            accVals[1][ivert] = (double)Vacc_ivdwAcc(thee, coord, radius);
        } else if (meth == 2) {
            accVals[0][ivert] = (double)Vacc_vdwAcc(thee, coord);
            accVals[1][ivert] = (double)Vacc_vdwAcc(thee, coord);
        } else VASSERT(0);
    }
    Gem_writeGMV(gm, iodev, iofmt, iohost, iofile, 1, accVals);
    Vmem_free(thee->vmem, Gem_numVV(gm), sizeof(double), 
      (void **)&(accVals[0]));
    Vmem_free(thee->vmem, Gem_numVV(gm), sizeof(double), 
      (void **)&(accVals[1]));
}
#endif /* defined(HAVE_FETK_H) */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_sphere
//
// Purpose:  Generates thee->npts somewhat uniformly distributed across a
//           sphere of unit radius centered at the origin.  Returns a
//           (npts x 3) double array, which you are resposible for destroying,
//           of approximatel the specified number of points; the actual number
//           is stored in the argument npts.  This routine was shamelessly
//           ripped off of sphere.F from UHBD as developed by Michael K. Gilson
//
// Author:   Nathan Baker (original FORTRAN routine from UHBD by Michael
//           Gilson)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double** Vacc_sphere(Vacc *thee, int *npts) {

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
    points = Vmem_malloc(thee->vmem, nactual, sizeof(double *));
    VASSERT(points != VNULL);
    for (i=0; i<nactual; i++) {
        points[i] = Vmem_malloc(thee->vmem, 3, sizeof(double));
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
                points[nactual][0] = cosphi * sintheta;
                points[nactual][1] = sinphi * sintheta;
                points[nactual][2] = costheta;
                nactual++;
            }
        }
    }

    *npts = nactual;
    return points;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_totalSASA
//
// Purpose:  Calculates the solvent-accessible area of the entire molecule
//
// Args:     radius  The radius of the solvent probe in Angstroms
//
// Author:   Nathan Baker (original FORTRAN routine from UHBD by Brock Luty)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vacc_totalSASA(Vacc *thee, double radius) { 

    int i;
    double area = 0.0;

    for (i=0; i<Valist_getNumberAtoms(thee->alist); i++) {
        thee->area[i] = Vacc_atomSASA(thee, radius, i);
        area += (thee->area[i]);
    }

    return area;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vacc_atomSASA
//
// Purpose:  Calculates the contribution to the PROBE-CENTERED
//           solvent-accessible area from this atom
//
// Args:     radius  The radius of the solvent probe in Angstroms
//           iatom   Index of the atom in thee->alist
//
// Author:   Nathan Baker (original FORTRAN routine from UHBD by Brock Luty)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vacc_atomSASA(Vacc *thee, double srad, int iatom) { 

    int ipt, covered;
    double area = 0.0;
    double *tPos, tRad, vec[3];
    Vatom *thisAtom;

    /* Get the atom information */
    thisAtom = Valist_getAtom(thee->alist, iatom);
    tPos = Vatom_getPosition(thisAtom);
    tRad = Vatom_getRadius(thisAtom);

    covered = 0;
    for (ipt=0; ipt<thee->nsphere; ipt++) {
        vec[0] = (tRad+srad)*thee->sphere[ipt][0] + tPos[0];
        vec[1] = (tRad+srad)*thee->sphere[ipt][1] + tPos[1];
        vec[2] = (tRad+srad)*thee->sphere[ipt][2] + tPos[2];
        if (ivdwAccExclus(thee, vec, srad, iatom)) area += 1.0;
    }

    /* We will return UHBD's asas2: probe-centered solvent-accessible surface
     * area */
    area = area/((double)(thee->nsphere))*4.0*VPI*(tRad+srad)*(tRad+srad);

    return area;

}
