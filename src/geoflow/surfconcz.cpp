#include "modules.h"
#include "Mat.h"
#include <valarray>
#include <vector>
#include <cmath>

#include <iostream>
#include <fstream>

using namespace std;


double left(const valarray<double>& pr, double ev){
    return floor( (pr - ev).min()/comdata.dcel )*comdata.dcel - ev;
}

double right(const valarray<double>& pr, double ev){
    return ceil( (pr + ev).max()/comdata.dcel )*comdata.dcel + ev;
}

void domainini(double xyzr[MAXATOMS][XYZRWIDTH], const int natm, const double extvalue){
    valarray<double> atom_x(natm), atom_y(natm), atom_z(natm), atom_r(natm);
    for(int i=0; i<natm; ++i){
        atom_x[i] = xyzr[i][0];
        atom_y[i] = xyzr[i][1];
        atom_z[i] = xyzr[i][2];
        atom_r[i] = xyzr[i][3];
    }
    
    double xleft = left(atom_x - atom_r, extvalue);
    double yleft = left(atom_y - atom_r, extvalue);
    double zleft = left(atom_z - atom_r, extvalue);

    double xright = right(atom_x + atom_r, extvalue);
    double yright = right(atom_y + atom_r, extvalue);
    double zright = right(atom_z + atom_r, extvalue);
    
    int nx = (xright - xleft)/comdata.dcel + 1;
    int ny = (yright - yleft)/comdata.dcel + 1;
    int nz = (zright - zleft)/comdata.dcel + 1;

    xright = xleft + comdata.dcel*(nx - 1);
    yright = yleft + comdata.dcel*(ny - 1);
    zright = zleft + comdata.dcel*(nz - 1);

    //set the stupid globals...
    comdata.xleft = xleft; comdata.xright = xright;
    comdata.yleft = yleft; comdata.yright = yright;
    comdata.zleft = zleft; comdata.zright = zright;
    comdata.nx = nx; comdata.ny = ny; comdata.nz = nz;
}

void volumintegration(const double* f, const int& nx, const int& ny, const int& nz, const double& dcel, double& volume){ //FIXME return volume
    valarray<double> ff(f, nx*ny*nz);//maybe divide before sum to avoid floating pt stuff
    ff /= 1000.0;
    ff *= dcel*dcel*dcel;
    volume = ff.sum();
}



double volumeIntegration(Mat<> f, double dcel){
    valarray<double> ff(f.data(), f.size());
    ff /= 1000.0;
    ff *= dcel*dcel*dcel;
    return ff.sum();
}

void upwinding(int nx, int ny, int nz, double dx, double dt, int nt, Mat<>& g, Mat<>& surfu, Mat<>& phitotx){
    double w1i[] = {-0.5, 0, 0.5};
    valarray<double> w1(w1i, 3);
    w1/=dx;
    
    double w2i[] = {1.0, -2.0, 1.0};
    valarray<double> w2(w2i, 3);
    w2/=dx*dx;
    
    double wxy = 0.25/dx/dx; //sign is sgn(prod(coord shifts))
    
    vector<double> surfnewi(surfu.data(), surfu.end()); //hide me in ctor
    Mat<> surfnew(surfnewi.data(), nx,ny,nz);
    
    for(int t=0; t<nt; ++t){ 
        for(int x=1; x<nx-1; ++x){ 
        for(int y=1; y<ny-1; ++y){ 
        for(int z=1; z<nz-1; ++z){
            if(g(x,y,z) > 2e-2){
                double phix = w1[0]*surfu(x-1,y,z) + w1[2]*surfu(x+1,y,z);
                double phiy = w1[0]*surfu(x,y-1,z) + w1[2]*surfu(x,y+1,z);
                double phiz = w1[0]*surfu(x,y,z-1) + w1[2]*surfu(x,y,z+1);
                
                double phixx = w2[0]*surfu(x-1,y,z) + w2[1]*surfu(x,y,z) + w2[2]*surfu(x+1,y,z);
                double phiyy = w2[0]*surfu(x,y-1,z) + w2[1]*surfu(x,y,z) + w2[2]*surfu(x,y+1,z);
                double phizz = w2[0]*surfu(x,y,z-1) + w2[1]*surfu(x,y,z) + w2[2]*surfu(x,y,z+1);
                
                double phixz = wxy*(surfu(x-1,y,z-1) + surfu(x+1,y,z+1)) - wxy*(surfu(x-1,y,z+1) + surfu(x+1,y,z-1));
                double phixy = wxy*(surfu(x-1,y-1,z) + surfu(x+1,y+1,z)) - wxy*(surfu(x-1,y+1,z) + surfu(x+1,y-1,z));
                double phiyz = wxy*(surfu(x,y-1,z-1) + surfu(x,y+1,z+1)) - wxy*(surfu(x,y-1,z+1) + surfu(x,y+1,z-1));
                
                double dphi =  (1.0 + phix*phix + phiy*phiy)*phizz
                             + (1.0 + phix*phix + phiz*phiz)*phiyy
                             + (1.0 + phiy*phiy + phiz*phiz)*phixx;
                
                dphi -= 2*( phix*phiy*phixy + phix*phiz*phixz + phiy*phiz*phiyz );            
                double gram = 1.0 + phix*phix + phiy*phiy + phiz*phiz;
                dphi = dphi/gram + sqrt(gram)*phitotx(x,y,z);
                
                surfnew(x,y,z) = min(1000.0, max(0.0, surfu(x,y,z) + dt*dphi));
            }
        }}}
        
        surfu = surfnew;
    }

}

void initial(double xl, double yl, double zl, double dx, int n_atom, const valarray<double>& atom_x, const valarray<double>& atom_y, const valarray<double>& atom_z, const valarray<double>& atom_r, Mat<>& g, Mat<>& phi){
    g = 1.0;
    phi = 0.0;

    double alpha = 1e3;
    
    for(int a=0; a<n_atom; ++a){
        double r = atom_r[a];
        double r2 = r*r;
        double zmin = ((atom_z[a] - zl - r)/dx + 1.0);
        double zmax = ((atom_z[a] - zl + r)/dx + 1.0);

        for(int z = ceil(zmin); z<=floor(zmax); ++z){
            double distxy = (zl + (z-1)*dx - atom_z[a]);
            double distxy2 = distxy*distxy;
            double rxy2 = abs(r2 - distxy2);
            double rxy = sqrt(rxy2);
            double ymin = ((atom_y[a] - yl - rxy)/dx + 1.0);
            double ymax = ((atom_y[a] - yl + rxy)/dx + 1.0);

            for(int y=ceil(ymin); y<=floor(ymax); ++y){
                double distx = (yl + (y-1)*dx - atom_y[a]);
                double distx2 = distx*distx;
                double rx = sqrt(abs(rxy2 - distx2));
                double xmin = ((atom_x[a] - xl - rx)/dx + 1.0);
                double xmax = ((atom_x[a] - xl + rx)/dx + 1.0);

                for(int x=ceil(xmin); x<=floor(xmax); ++x){
                    g(x,y,z) = 0;
                    phi(x,y,z) = alpha;
                }
            }
        }
    }
}

void potIntegral(double rcfactor, double ddx, int natm, valarray<double>& atom_x,valarray<double>& atom_y,valarray<double>& atom_z,valarray<double>& seta12, valarray<double>& seta6, valarray<double>& epsilon, valarray<double>& sigma, Mat<>& g, Mat<>& potr, Mat<>& pota){
    for(int x=2; x<potr._nx; ++x){ 
    for(int y=2; y<potr._ny; ++y){ 
    for(int z=2; z<potr._nz; ++z){
        if(g(x,y,z) == 0){ continue; }
        double pr=0, pa=0;
        for(int a=0; a<natm; ++a){
            double xi = comdata.xleft + (x-1)*ddx;
            double yi = comdata.yleft + (y-1)*ddx;
            double zi = comdata.zleft + (z-1)*ddx;
            double dist = sqrt( dot(xi-atom_x[a], yi-atom_y[a], zi-atom_z[a]) ) + lj.prob;
            double ratio = (dist==0.0) ? 1.0 : sigma[a]/dist;
    
            if(lj.iwca == 1){
                if(ratio*rcfactor > 1){
                    pr += seta12[a]*pow(ratio,12.0) - seta6[a]*pow(ratio,6.0) + epsilon[a];
                    pa -= epsilon[a];
                }else{
                    pa = pa - seta6[a]*pow(ratio,6.0) + seta12[a]*pow(ratio,12.0);
                }
            }else{
                pr += pow(ratio,12.0)*seta12[a];
                pa -= pow(ratio, 6.0)*seta6[a];
            }
            potr(x,y,z) = pr; 
            pota(x,y,z) = pa;
        }
    }}}
}

//void adicor(double* su, double* g, double& ddx, double& dt, int& nz,int& ny,int& nx, double* phitotx);

void yhsurface(double xyzr[MAXATOMS][XYZRWIDTH], double* ljepsilon, int natm, double tott, double dt, double* _phitotx, double* _surfu, int iloop, double& area, double& volume, double& attint, double alpha, int iadi, int igfin
){
    int nx = comdata.nx, ny = comdata.ny, nz = comdata.nz;
    double xl = comdata.xleft, yl = comdata.yleft, zl = comdata.zleft;
    double ddx = comdata.dcel;
    valarray<double> atom_x(natm), atom_y(natm), atom_z(natm), atom_r(natm);
    for(int i=0; i<natm; ++i){
        atom_x[i] = xyzr[i][0];
        atom_y[i] = xyzr[i][1];
        atom_z[i] = xyzr[i][2];
        atom_r[i] = xyzr[i][3];
    }

    vector<double> _su(nx*ny*nz);
    Mat<> su(_su, nx,ny,nz);
    vector<double> _g(nx*ny*nz);
    Mat<> g(_g, nx,ny,nz);

    initial(xl,yl,zl, ddx, natm, atom_x,atom_y,atom_z,atom_r, g, su);

    Mat<> surfu(_surfu, nx,ny,nz);
    if(iloop > 1 && igfin == 1){ su = surfu; }


    valarray<double> sigma(natm), seta12(natm), seta6(natm), expan(natm), epsilon(natm);
    double rcfactor = (lj.ffmodel == 1) ? 1.0 : pow(2.0, 1.0/6.0);
    if(lj.ffmodel == 1){
        for(int i=0; i<natm; ++i){
            sigma[i] = atom_r[i] + lj.sigmas;
            expan[i] = atom_r[i] + lj.prob; //doesnt need to be a vector
            if(lj.vdwdispersion == 0){
                epsilon[i] = 0;
                seta12[i] = 0;
                seta6[i] = 0;
            }else{
                double se = sigma[i]/expan[i];
                epsilon[i] = pow( pow(se, 12.0) - 2.0*pow(se, 6.0) , -1.0);
                seta12[i] = lj.iosetar*lj.vdwdispersion*epsilon[i];
                seta6[i] = 2.0*lj.iosetaa*lj.vdwdispersion*epsilon[i];
            }
        }
    }else{
        for(int i=0; i<natm; ++i){
            sigma[i] = sqrt(4.0*atom_r[i]*lj.sigmas);
            if(lj.vdwdispersion == 0){
                epsilon[i] = 0;
                seta12[i] = 0;
                seta6[i] = 0;
            }else{
                epsilon[i] = sqrt(ljepsilon[i]*lj.epsilonw);
                seta12[i] = 4.0*epsilon[i];
                seta6[i] = 4.0*epsilon[i];
            }
        }
    }

    vector<double> _potr(nx*ny*nz), _pota(nx*ny*nz);
    Mat<> potr(_potr, nx,ny,nz), pota(_pota, nx,ny,nz);
    potr=0.0, pota=0.0;
    potIntegral(rcfactor, ddx, natm, atom_x, atom_y, atom_z, seta12, seta6, epsilon, sigma, g, potr, pota);

    if(lj.iwca == 1) potr = 0;

    Mat<> phitotx(_phitotx, nx,ny,nz);
    for(int i=0; i<phitotx.size(); ++i){ 
        phitotx[i] = -lj.conms - phitotx[i] + lj.roro*(potr[i] + pota[i]);
    }

    if(iadi==0 || iloop>1){
        int nt = ceil(tott/dt);
        dt = ddx*ddx/4.5;
        upwinding(nx,ny,nz, ddx,dt, nt, g,su,phitotx);
    }else{
        cerr << "ADI not implemented..." << endl;
        exit(1);
        //for(int t=1; t<=nt; ++t){  dt is unset in parent call!!!
        //    adicor(su.data(), g.data(), ddx, dt, nz,ny,nx, phitotx.data()); //this has issues, should have &
        //}
    }

    if(iloop < 2){
        surfu = su;
    }else{
        for(int i=0; i<surfu.size(); ++i){ 
           surfu[i] = surfu[i]*alpha + su[i]*(1.0 - alpha);
        }
        su = surfu;
    }

    volume = volumeIntegration(su, ddx);

    vector<double> _fintegr(nx*ny*nz);
    Mat<> fintegr(_fintegr, nx,ny,nz);
    double weight = pow(2.0*ddx, -1.0);
    for(int x=2; x<nx; ++x){
    for(int y=2; y<ny; ++y){
    for(int z=2; z<nz; ++z){
        double sux = su(x+1,y,z) - su(x-1,y,z);
        double suy = su(x,y+1,z) - su(x,y-1,z);
        double suz = su(x,y,z+1) - su(x,y,z-1);
        fintegr(x,y,z) = sqrt(dot(sux,suy,suz))*weight;
    }}}
    area = volumeIntegration(fintegr, ddx);
   
    potIntegral(rcfactor, ddx, natm, atom_x, atom_y, atom_z, seta12, seta6, epsilon, sigma, g, potr, pota);

    if(lj.iwca==1){
        for(int i=0; i<fintegr.size(); ++i){
            fintegr[i] = pota[i]*(1e3 - su[i]);
        }
    }else{
        for(int i=0; i<fintegr.size(); ++i){ 
            fintegr[i] = (pota[i] + potr[i])*(1e3 - su[i]);
        }

    }

    attint = volumeIntegration(fintegr, ddx);
}

