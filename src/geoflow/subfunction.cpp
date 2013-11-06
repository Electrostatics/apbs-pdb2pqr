#include "modules.h"
#include "Mat.h"
#include <valarray>
#include <cmath>

using namespace std;

extern "C"{

double qbboundary(int& natm, double& x,double& y,double& z, double xyzr[MAXATOMS][XYZRWIDTH], double* pqr, double& epsilonsp){
    double vbdn = 0;
    for(int a=0; a<natm; ++a){
        double x_q = x - xyzr[a][1];
        double y_q = y - xyzr[a][2];
        double z_q = z - xyzr[a][3];
        double q_q = pqr[a];
        double rr = sqrt(dot(x_q,y_q,z_q));
        vbdn += q_q/(epsilonsp*rr);
    }
    return vbdn;
}


double qbinterior(int& natm, double& x,double& y,double& z, double* _charget, double* _corlocqt){
    Mat<> charget(_charget, natm,8), corlocqt(_corlocqt, natm,8,3);
    double fp = 0;
    for(int a=1; a<=natm; ++a){
    for(int ii=1; ii<=8; ++ii){
        double xc = x - corlocqt(a,ii,1);
        double yc = y - corlocqt(a,ii,2);
        double zc = z - corlocqt(a,ii,3);
        if( dot(xc,yc,zc) <= 1e-13){
            fp -= 4.0*comdata.pi*charget(a,ii)/(comdata.deltax*comdata.deltay*comdata.deltaz);
        }
    }}
    return fp;
}

double qb(int& natm,int& i,int& j,int& k, double xyzr[MAXATOMS][XYZRWIDTH], double*pqr, double* charget, double* corlocqt, double& epsilonsp){
    double x = xvalue(i);
    double y = yvalue(j);
    double z = zvalue(k);
    if(i<2 || i>comdata.nx-1 || j<2 || j>comdata.ny-1 || k<2 || k>comdata.nz-1){
        return qbboundary(natm, x,y,z, xyzr, pqr, epsilonsp);
    }else{
        return qbinterior(natm,x,y,z, charget, corlocqt);
    }
}

void seteqb(double* bg, double xyzr[MAXATOMS][XYZRWIDTH], double* pqr, int& natm, double* charget, double* corlocqt, double* epsilonsp){
    for(int i=1; i<=comdata.nx; ++i){
    for(int j=1; j<=comdata.ny; ++j){
    for(int k=1; k<=comdata.nz; ++k){
        double fp = qb(natm,i,j,k,xyzr,pqr,charget,corlocqt,*epsilonsp);
        int ijk = (i-1)*comdata.ny*comdata.nz + (j-1)*comdata.nz + k - 1;
        bg[ijk] = fp;
        
    }}}
}

}//extern
