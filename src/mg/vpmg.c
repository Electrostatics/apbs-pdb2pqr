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
// File:     vpmg.c
//
// Purpose:  Class Vpmg: methods.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "apbs/vpmg.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmg: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VACC)
#endif /* if !defined(VINLINE_VACC) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmg: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmg: Private methods
/////////////////////////////////////////////////////////////////////////// */
#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define IJKx(j,k,i) (((i)*(ny)*(nz))+((k)*(ny))+(j))
#define IJKy(i,k,j) (((j)*(nx)*(nz))+((k)*(nx))+(i))
#define IJKz(i,j,k) (((k)*(nx)*(ny))+((j)*(nx))+(i))


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  MGpde_anaU
//
// Purpose:  Dirichlet boundary function and initial approximation function.
//
// Args:     x    = position vector
//           flag = evaluation flag 
//                    0 => single Debye-Huckel sphere
//                    1 => Debye-Huckel sphere for each atom
//
// Author:   Nathan Baker and Michael Holst
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double myGreens(Vpbe *pbe, double x[], int flag) {

    double size, *position, charge, xkappa, eps_w, dist, T, val, pot;
    int i, iatom;
    Vatom *atom;
    Valist *alist;

    if (flag == 0) {
        /* Get the solute radius in meters and position in angstroms */
        size = (1.0e-10)*Vpbe_getSoluteRadius(pbe);
        position = Vpbe_getSoluteCenter(pbe);

        /* We keep the charge relative to units of ec that are factored out;
         * this term should be dimensionless. The dielectric is unitless. */
        charge = Vunit_ec*Vpbe_getSoluteCharge(pbe);
        eps_w = Vpbe_getSolventDiel(pbe);

        /* Get xkappa in units of inverse meters */
        xkappa = (1.0e10)*Vpbe_getXkappa(pbe);

        /* The temperature is in units of K */
        T = Vpbe_getTemperature(pbe);

        /* Compute the distance (in units of m) */
        dist = 0;
        for (i=0; i<3; i++) 
            dist += ((position[i] - x[i])*(position[i] - x[i])); 
        dist = (1.0e-10)*VSQRT(dist);

        /* Compute the potential in J/electron */
        val = (charge)/(4*VPI*Vunit_eps0*eps_w*dist);
        if (xkappa != 0.0) val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
        /* Scale the potential to be dimensionless */
        val = val*Vunit_ec/(Vunit_kb*T);
        return val;
    } else if (flag == 1) {
        pot = 0.0;
        eps_w = Vpbe_getSolventDiel(pbe);
        xkappa = (1.0e10)*Vpbe_getXkappa(pbe);
        T = Vpbe_getTemperature(pbe);
        alist = Vpbe_getValist(pbe);
        for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
            atom = Valist_getAtom(alist, iatom);
            position = Vatom_getPosition(atom);
            charge = Vunit_ec*Vatom_getCharge(atom);
            size = (1e-10)*Vatom_getRadius(atom);
            dist = 0;
            for (i=0; i<3; i++)
               dist += ((position[i] - x[i])*(position[i] - x[i]));
            dist = (1.0e-10)*VSQRT(dist);
            val = (charge)/(4*VPI*Vunit_eps0*eps_w*dist);
            if (xkappa != 0.0)
              val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
            val = val*Vunit_ec/(Vunit_kb*T);
            pot = pot + val;
        }

        return pot;
    } 

    return 0;
}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_ctor
//
// Purpose:  Construct the PMG parameter object; see header file for more
//           information
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vpmg* Vpmg_ctor(Vpmgp *pmgp, Vpbe *pbe) {

    Vpmg *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vpmg) );
    VASSERT( thee != VNULL);
    VASSERT(Vpmg_ctor2(thee, pmgp, pbe));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_ctor2
//
// Purpose:  Construct the PMG parameter object
//
// Notes:    See header files for default parameter values
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpmg_ctor2(Vpmg *thee, Vpmgp *pmgp, Vpbe *pbe) {

    int nxc, nyc, nzc, nf, nc, narrc, n_rpc, n_iz, n_ipc;

    /* Get the parameters */    
    VASSERT(pmgp != VNULL); 
    VASSERT(pbe != VNULL); 
    thee->pmgp = pmgp;
    thee->pbe = pbe;

    if (thee->pmgp->nonlin != 0) {
        Vnm_print(2, "Vpmg_ctor2: Sorry, nonlinear PBE support not available yet.\n");
        VASSERT(0);
    }

    /* Set up the memory */
    thee->vmem = Vmem_ctor("APBS:VPMG");

    /* Calculate storage requirements */
    mgsz_(&(thee->pmgp->mgcoar), &(thee->pmgp->mgdisc), &(thee->pmgp->mgsolv),
      &(thee->pmgp->nx), &(thee->pmgp->ny), &(thee->pmgp->nz),
      &(thee->pmgp->nlev), &nxc, &nyc, &nzc, &nf, &nc, &(thee->pmgp->narr),
      &narrc, &n_rpc, &n_iz, &n_ipc, &(thee->pmgp->nrwk), &(thee->pmgp->niwk));

    /* Allocate storage */
    thee->iparm = (int *)Vmem_malloc(thee->vmem, 100, sizeof(int));
    thee->rparm = (double *)Vmem_malloc(thee->vmem, 100, sizeof(double));
    thee->iwork = (int *)Vmem_malloc(thee->vmem, thee->pmgp->niwk, 
      sizeof(int));
    thee->rwork = (double *)Vmem_malloc(thee->vmem, thee->pmgp->nrwk, 
      sizeof(double));
    thee->a1cf = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr, 
      sizeof(double));
    thee->a2cf = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr, 
      sizeof(double));
    thee->a3cf = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr, 
      sizeof(double));
    thee->ccf = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr, 
      sizeof(double));
    thee->fcf = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr, 
      sizeof(double));
    thee->tcf = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr, 
      sizeof(double));
    thee->u = (double *)Vmem_malloc(thee->vmem, thee->pmgp->narr, 
      sizeof(double));
    thee->xf = (double *)Vmem_malloc(thee->vmem, 5*(thee->pmgp->nx), 
      sizeof(double));
    thee->yf = (double *)Vmem_malloc(thee->vmem, 5*(thee->pmgp->ny), 
      sizeof(double));
    thee->zf = (double *)Vmem_malloc(thee->vmem, 5*(thee->pmgp->nz), 
      sizeof(double));
    thee->gxcf = (double *)Vmem_malloc(thee->vmem, 
      10*(thee->pmgp->ny)*(thee->pmgp->nz), sizeof(double));
    thee->gycf = (double *)Vmem_malloc(thee->vmem, 
      10*(thee->pmgp->nx)*(thee->pmgp->nz), sizeof(double));
    thee->gzcf = (double *)Vmem_malloc(thee->vmem, 
      10*(thee->pmgp->nx)*(thee->pmgp->ny), sizeof(double));

    /* Plop some of the parameters into the iparm and rparm arrays */
    packmg_(thee->iparm, thee->rparm, &(thee->pmgp->nrwk), &(thee->pmgp->niwk),
      &(thee->pmgp->nx), &(thee->pmgp->ny), &(thee->pmgp->nz),
      &(thee->pmgp->nlev), &(thee->pmgp->nu1), &(thee->pmgp->nu2),
      &(thee->pmgp->mgkey), &(thee->pmgp->itmax), &(thee->pmgp->istop),
      &(thee->pmgp->ipcon), &(thee->pmgp->nonlin), &(thee->pmgp->mgsmoo),
      &(thee->pmgp->mgprol), &(thee->pmgp->mgcoar), &(thee->pmgp->mgsolv),
      &(thee->pmgp->mgdisc), &(thee->pmgp->iinfo), &(thee->pmgp->errtol),
      &(thee->pmgp->ipkey), &(thee->pmgp->omegal), &(thee->pmgp->omegan),
      &(thee->pmgp->irite), &(thee->pmgp->iperf));

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_solve
//
// Purpose:  Solve the equation
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_solve(Vpmg *thee) {

    switch(thee->pmgp->meth) {
        /* CGMG (linear) */
        case 0:
            cgmgdriv_(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
              thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
              thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            break;
        /* Newton (nonlinear) */
        case 1:
            newdriv_(thee->iparm, thee->rparm, thee->iwork, thee->rwork, 
              thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
              thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf, 
              thee->fcf, thee->tcf);
            break;
        /* MG (linear/nonlinear) */
        case 2:
	    mgdriv_(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            break;
        /* CGHS (linear/nonlinear) */
        case 3: 
	    ncghsdriv_(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            break;
        /* SOR (linear/nonlinear) */
        case 4:
	    nsordriv_(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            break;
        /* GSRB (linear/nonlinear) */
        case 5:
	    ngsrbdriv_(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf); 
            break;
        /* WJAC (linear/nonlinear) */
        case 6:
	    nwjacdriv_(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            break;
        /* RICH (linear/nonlinear) */
        case 7:
	    nrichdriv_(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
	      thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
	      thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
              thee->fcf, thee->tcf);
            break;
        /* Error handling */
        default: 
            Vnm_print(2, "Vpgm_solve: invalid solver method key (%d)\n",
              thee->pmgp->key);
            break;
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_fillco
//
// Purpose:  Fill the coefficient arrays prior to solving the equation
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_fillco(Vpmg *thee) {

    Vacc *acc;
    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xlen, ylen, zlen, position[3], ifloat, jfloat, kfloat, accf;
    double zmagic, irad, srad, charge, dx, dy, dz, zkappa2, epsw, epsp;
    double hx, hy, hz;
    int i, j, k, nx, ny, nz, iatom, ihi, ilo, jhi, jlo, khi, klo;
    int acclo, accmid, acchi;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    irad = Vpbe_getIonRadius(pbe);
    srad = Vpbe_getSolventRadius(pbe);
    zmagic = Vpbe_getZmagic(pbe);
    zkappa2 = Vpbe_getZkappa2(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    epsp = Vpbe_getSoluteDiel(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hz = thee->pmgp->hz;
   
    /* Define the total domain size */
    xlen = hx*(nx - 1);
    ylen = hy*(ny - 1);
    zlen = hz*(nz - 1);

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);
    thee->rparm[2] = xmin;
    thee->rparm[3] = xmax;
    thee->rparm[4] = ymin;
    thee->rparm[5] = ymax;
    thee->rparm[6] = zmin;
    thee->rparm[7] = zmax;

    /* Fill the mesh point coordinate arrays */
    for (i=0; i<nx; i++) thee->xf[i] = xmin + i*hx;
    for (i=0; i<ny; i++) thee->yf[i] = ymin + i*hy;
    for (i=0; i<nz; i++) thee->zf[i] = zmin + i*hz;

    /* Fill the coefficient arrays */
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {

                position[0] = thee->xf[i];
                position[1] = thee->yf[j];
                position[2] = thee->zf[k];

                /* the scalar (0th derivative) entry */
                if (Vacc_ivdwAcc(acc, position, irad) == 1) 
                  thee->ccf[IJK(i,j,k)] = zkappa2;
                else thee->ccf[IJK(i,j,k)] = 0.0;

                /* the chosen true solution */
                thee->tcf[IJK(i,j,k)] = 0.0;

                /* Clear out the load vector */
                thee->fcf[IJK(i,j,k)] = 0.0;

                /* The diagonal tensor (2nd derivative) entries.  Each of these
                 * entries is evaluated ad the grid edges midpoints.  We will
                 * implement a very rudimentary form of dielectric smoothing.
                 * Specifically, the dielectric will be evaluated at the mid
                 * point and the two flanking mesh points.  The fraction of the
                 * grid edge in the solvent will then be calculated from these
                 * three values (i.e., either 0, 1/3, 2/3, or 1).  The
                 * dielectric value at the midpoint will then be assigned based
                 * on the usual dielectric smoothing formula:
                 * \epsilon_s\epsilon_i/(a\epsilon_s + (1-a)\epsilon_i)  */
                /* x-direction */
                position[0] = thee->xf[i] + hx/2.0;
                position[1] = thee->yf[j];
                position[2] = thee->zf[k];
                accmid = Vacc_molAcc(acc, position, srad);
                position[0] = thee->xf[i];
                acclo = Vacc_molAcc(acc, position, srad);
                position[0] = thee->xf[i] + hx;
                acchi = Vacc_molAcc(acc, position, srad);
                accf = ((double)acchi + (double)accmid + (double)acclo)/3.0;
                thee->a1cf[IJK(i,j,k)] = 
                  epsw*epsp/((1-accf)*epsw + accf*epsp);
                /* y-direction */
                position[0] = thee->xf[i];
                position[1] = thee->yf[j] + hy/2.0;
                position[2] = thee->zf[k];
                accmid = Vacc_molAcc(acc, position, srad);
                position[1] = thee->yf[j];
                acclo = Vacc_molAcc(acc, position, srad);
                position[1] = thee->yf[j] + hy;
                acchi = Vacc_molAcc(acc, position, srad);
                accf = ((double)acchi + (double)accmid + (double)acclo)/3.0;
                thee->a2cf[IJK(i,j,k)] = 
                  epsw*epsp/((1-accf)*epsw + accf*epsp);
                /* y-direction */
                position[0] = thee->xf[i];
                position[1] = thee->yf[j];
                position[2] = thee->zf[k] + hz/2.0;
                accmid = Vacc_molAcc(acc, position, srad);
                position[2] = thee->zf[k];
                acclo = Vacc_molAcc(acc, position, srad);
                position[2] = thee->zf[k] + hz;
                acchi = Vacc_molAcc(acc, position, srad);
                accf = ((double)acchi + (double)accmid + (double)acclo)/3.0;
                thee->a3cf[IJK(i,j,k)] = 
                  epsw*epsp/((1-accf)*epsw + accf*epsp);
            }
        }
    }

    /* Fill source term array */
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
        atom = Valist_getAtom(alist, iatom);

        position[0] = Vatom_getPosition(atom)[0];
        position[1] = Vatom_getPosition(atom)[1];
        position[2] = Vatom_getPosition(atom)[2];

        /* Make sure we're on the grid */
        if ((position[0]<=xmin) || (position[0]>=xmax)  || \
            (position[1]<=ymin) || (position[1]>=ymax)  || \
            (position[2]<=zmin) || (position[2]>=zmax)) {
            Vnm_print(2, "MGpde_fillco:  Atom #%d at (%4.3f, %4.3f, %4.3f) is
off the mesh (ignoring)!\n",
            iatom, position[0], position[1], position[2]);
        } else {

            /* Figure out which vertices we're next to */
            ifloat = (position[0] - xmin)/(hx);
            jfloat = (position[1] - xmin)/(hy);
            kfloat = (position[2] - xmin)/(hz);

            ihi = (int)ceil(ifloat);
            ilo = (int)floor(ifloat);
            jhi = (int)ceil(jfloat);
            jlo = (int)floor(jfloat);
            khi = (int)ceil(kfloat);
            klo = (int)floor(kfloat);

            /* Now assign fractions of the charge to the nearby verts */
            charge = zmagic*Vatom_getCharge(atom)/hx/hy/hz;
            dx = ifloat - (double)(ilo);
            dy = jfloat - (double)(jlo);
            dz = kfloat - (double)(klo);
            thee->fcf[IJK(ihi,jhi,khi)] += (dx*dy*dz*charge);
            thee->fcf[IJK(ihi,jlo,khi)] += (dx*(1.0-dy)*dz*charge);
            thee->fcf[IJK(ihi,jhi,klo)] += (dx*dy*(1.0-dz)*charge);
            thee->fcf[IJK(ihi,jlo,klo)] += (dx*(1.0-dy)*(1.0-dz)*charge);
            thee->fcf[IJK(ilo,jhi,khi)] += ((1.0-dx)*dy*dz*charge);
            thee->fcf[IJK(ilo,jlo,khi)] += ((1.0-dx)*(1.0-dy)*dz*charge);
            thee->fcf[IJK(ilo,jhi,klo)] += ((1.0-dx)*dy*(1.0-dz)*charge);
            thee->fcf[IJK(ilo,jlo,klo)] += ((1.0-dx)*(1.0-dy)*(1.0-dz)*charge);

        }
    }

    /* the "i" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            position[0] = thee->xf[0];
            position[1] = thee->yf[j];
            position[2] = thee->zf[k];
            thee->gxcf[IJKx(j,k,0)] = myGreens(pbe, position, 
              thee->pmgp->bcfl);
            position[0] = thee->xf[nx-1];
            thee->gxcf[IJKx(j,k,1)] = myGreens(pbe, position, 
              thee->pmgp->bcfl);
            thee->gxcf[IJKx(j,k,2)] = 0.0;
            thee->gxcf[IJKx(j,k,3)] = 0.0;
        }
    }

    /* the "j" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        for (i=0; i<nx; i++) {
            position[0] = thee->xf[i];
            position[1] = thee->yf[0];
            position[2] = thee->zf[k];
            thee->gycf[IJKy(i,k,0)] = myGreens(pbe, position, 
              thee->pmgp->bcfl);
            position[1] = thee->yf[ny-1];
            thee->gycf[IJKy(i,k,1)] = myGreens(pbe, position, 
              thee->pmgp->bcfl);
            thee->gycf[IJKy(i,k,2)] = 0.0;
            thee->gycf[IJKy(i,k,3)] = 0.0;
        }
    }

    /* the "k" boundaries (dirichlet) */
    for (j=0; j<ny; j++) {
        for (i=0; i<nx; i++) {
            position[0] = thee->xf[i];
            position[1] = thee->yf[j];
            position[2] = thee->zf[0];
            thee->gzcf[IJKz(i,j,0)] = myGreens(pbe, position, 
              thee->pmgp->bcfl);
            position[2] = thee->zf[nz-1];
            thee->gzcf[IJKz(i,j,1)] = myGreens(pbe, position, 
              thee->pmgp->bcfl);
            thee->gzcf[IJKz(i,j,2)] = 0.0;
            thee->gzcf[IJKz(i,j,3)] = 0.0;
        }
    }
}
    
/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_getLinearEnergy1
//
// Purpose:  using the solution at the finest mesh level, get the
//           electrostatic energy using the free energy functional for the
//           linearized Poisson-Boltzmann equation without removing any
//           self-interaction terms (i.e., removing the reference state of
//           isolated charges present in an infinite dielectric continuum with
//           the same relative permittivity as the interior of the protein).
//           In other words, we calculate
//             \[ G = \frac{1}{2} \sum_i q_i u(r_i) \]
//           and return the result in units of $k_B T$.  
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpmg_getLinearEnergy1(Vpmg *thee) {

    int iatom, nx, ny, nz, ihi, ilo, jhi, jlo, khi, klo;
    double xmax, xmin, ymax, ymin, zmax, zmin, hx, hy, hz, ifloat, jfloat;
    double charge, kfloat, dx, dy, dz, energy, uval, *position;
    Valist *alist;
    Vatom *atom; 
    Vpbe *pbe;

    pbe = thee->pbe;
    alist = pbe->alist;
    VASSERT(alist != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hz = thee->pmgp->hz;
    xmax = thee->xf[nx-1];
    ymax = thee->yf[ny-1];
    zmax = thee->zf[nz-1];
    xmin = thee->xf[0];
    ymin = thee->yf[0];
    zmin = thee->zf[0];
  
    energy = 0.0;

    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
        /* Get atomic information */
        atom = Valist_getAtom(alist, iatom);
        charge = Vatom_getCharge(atom);
        position = Vatom_getPosition(atom);

        /* Make sure the atom is on the grid */
        if ((position[0]<=xmin) || (position[0]>=xmax)  || \
            (position[1]<=ymin) || (position[1]>=ymax)  || \
            (position[2]<=zmin) || (position[2]>=zmax)) {
            Vnm_print(2, "MGpde_fillco:  Atom #%d at (%4.3f, %4.3f, %4.3f)
is off the mesh (ignoring)!\n",
                iatom, position[0], position[1], position[2]);
        } else {
            /* Figure out which vertices we're next to */
            ifloat = (position[0] - xmin)/hx;
            jfloat = (position[1] - xmin)/hy;
            kfloat = (position[2] - xmin)/hz;
            ihi = (int)ceil(ifloat);
            ilo = (int)floor(ifloat);
            jhi = (int)ceil(jfloat);
            jlo = (int)floor(jfloat);
            khi = (int)ceil(kfloat);
            klo = (int)floor(kfloat);

            /* Now get trilinear interpolation constants */
            dx = ifloat - (double)(ilo);
            dy = jfloat - (double)(jlo);
            dz = kfloat - (double)(klo);
            uval =  dx*dy*dz*(thee->u[IJK(ihi,jhi,khi)])
                  + dx*(1.0-dy)*dz*(thee->u[IJK(ihi,jlo,khi)])
                  + dx*dy*(1.0-dz)*(thee->u[IJK(ihi,jhi,klo)])
                  + dx*(1.0-dy)*(1.0-dz)*(thee->u[IJK(ihi,jlo,klo)])
                  + (1.0-dx)*dy*dz*(thee->u[IJK(ilo,jhi,khi)])
                  + (1.0-dx)*(1.0-dy)*dz*(thee->u[IJK(ilo,jlo,khi)])
                  + (1.0-dx)*dy*(1.0-dz)*(thee->u[IJK(ilo,jhi,klo)])
                  + (1.0-dx)*(1.0-dy)*(1.0-dz)*(thee->u[IJK(ilo,jlo,klo)]);
            energy += (uval*charge);
        }
    }

    energy = 0.5*energy;
    return energy;
}
    
/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_dtor
//
// Purpose:  Clean up
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_dtor(Vpmg **thee) {
    
    if ((*thee) != VNULL) {
        Vpmg_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vpmg), (void **)thee);
        (*thee) = VNULL;
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmg_dtor2
//
// Purpose:  Clean up
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmg_dtor2(Vpmg *thee) { 

    /* Clean up the storage */
    Vmem_free(thee->vmem, 100, sizeof(int), (void **)&(thee->iparm));
    Vmem_free(thee->vmem, 100, sizeof(double), (void **)&(thee->rparm));
    Vmem_free(thee->vmem, thee->pmgp->niwk, sizeof(int), 
      (void **)&(thee->iwork));
    Vmem_free(thee->vmem, thee->pmgp->nrwk, sizeof(double), 
      (void **)&(thee->rwork));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->a1cf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double), 
      (void **)&(thee->a2cf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->a3cf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->ccf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double), 
      (void **)&(thee->fcf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double), 
      (void **)&(thee->tcf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double), 
      (void **)&(thee->u));
    Vmem_free(thee->vmem, 5*(thee->pmgp->nx), sizeof(double),
      (void **)&(thee->xf));
    Vmem_free(thee->vmem, 5*(thee->pmgp->ny), sizeof(double),
      (void **)&(thee->yf));
    Vmem_free(thee->vmem, 5*(thee->pmgp->nz), sizeof(double),
      (void **)&(thee->zf));
    Vmem_free(thee->vmem, 10*(thee->pmgp->ny)*(thee->pmgp->nz), sizeof(double),
      (void **)&(thee->gxcf));
    Vmem_free(thee->vmem, 10*(thee->pmgp->nx)*(thee->pmgp->nz), sizeof(double),
      (void **)&(thee->gycf));
    Vmem_free(thee->vmem, 10*(thee->pmgp->nx)*(thee->pmgp->ny), sizeof(double),
      (void **)&(thee->gzcf));

    Vmem_dtor(&(thee->vmem));
}


