#include "modules.h"
#include "Mat.h"
#include <valarray>
#include <cmath>

using namespace std;

extern "C"{

int inverx(double& x){ return int( (x - comdata.xleft)/comdata.deltax ) + 1; }
int invery(double& y){ return int( (y - comdata.yleft)/comdata.deltay ) + 1; }
int inverz(double& z){ return int( (z - comdata.zleft)/comdata.deltaz ) + 1; }

void chargedist(double* _atmpos, double* chratm, int& natm, double* _charget, double* _corlocqt, int* _loc_qt, int& iatm){
    Mat<> atmpos(_atmpos, 4,natm);
    double x_q = atmpos(1, iatm);
    double y_q = atmpos(2, iatm);
    double z_q = atmpos(3, iatm);
    double q_q = chratm[iatm-1];

    int i_q = inverx(x_q);
    int j_q = invery(y_q);
    int k_q = inverz(z_q);


    vector<int> _loc_q(8*3);
    Mat<int> loc_q(_loc_q, 8,3);
    vector<double> _corlocq(8*3);
    Mat<> corlocq(_corlocq, 8,3);
    for(int i=0; i<=1; ++i){
    for(int j=0; j<=1; ++j){
    for(int k=0; k<=1; ++k){
        int ind = 4*k + 2*j + i + 1;
    
        int ind_i = i_q + i;
        int ind_j = j_q + j;
        int ind_k = k_q + k;
        
        corlocq(ind,1) = xvalue(ind_i);
        corlocq(ind,2) = yvalue(ind_j);
        corlocq(ind,3) = zvalue(ind_k);
        
        loc_q(ind,1) = ind_i;
        loc_q(ind,2) = ind_j;
        loc_q(ind,3) = ind_k;

    }}}

    double x = xvalue(i_q);
    double y = yvalue(j_q);
    double z = zvalue(k_q);

    double xd1 = x_q - x;
    double yd1 = y_q - y;
    double zd1 = z_q - z;

    valarray<int> iloc(1.0, 3);//ones
    if(xd1==0){ iloc[1]=0; }
    if(yd1==0){ iloc[2]=0; }
    if(zd1==0){ iloc[3]=0; }

    valarray<double> charge(0.0, 8);
    switch(iloc.sum()){
        case 3:
            for(int i=0; i<=1; ++i){
            for(int j=0; j<=1; ++j){
            for(int k=0; k<=1; ++k){
                int ind = 4*k + 2*j + i;
            
                double xd = i*comdata.deltax - xd1;
                double yd = j*comdata.deltay - yd1;
                double zd = k*comdata.deltaz - zd1;

                charge[ind] = 1.0/abs(xd*yd*zd);
            }}}
            break;

        case 2:
            if(iloc[0] == 0){
                for(int j=0; j<=1; ++j){
                for(int k=0; k<=1; ++k){
                    double yd = j*comdata.deltay - yd1;
                    double zd = k*comdata.deltaz - zd1;
                    charge[4*j + 2*k] = 1.0/abs(yd*zd);
                }}
            }else if(iloc[1] == 0){
                for(int i=0; i<=1; ++i){
                for(int k=0; k<=1; ++k){
                    double xd = i*comdata.deltax - xd1;
                    double zd = k*comdata.deltaz - zd1;
                    charge[i + 4*k] = 1.0/abs(xd*zd);
                }}
            }else if(iloc[2] == 0){
                for(int i=0; i<=1; ++i){
                for(int j=0; j<=1; ++j){
                    double xd = i*comdata.deltax - xd1;
                    double yd = j*comdata.deltay - yd1;
                    charge[i + 2*j] = 1.0/abs(xd*yd);
                }}
            }
            break;

        case 1:
            if(iloc[0]){
                charge[0] = 1.0/xd1;
                charge[1] = 1.0/(comdata.deltax-xd1);
            }else if(iloc[1]){
                charge[0] = 1.0/yd1;
                charge[2] = 1.0/(comdata.deltay-yd1);
            }else if(iloc[2]){
                charge[0] = 1.0/zd1;
                charge[4] = 1.0/(comdata.deltaz-zd1);
            }
            break;
       default:
           charge[0] = 1.0;
    }

    charge=q_q*charge/charge.sum();

    Mat<> charget(_charget, natm,8), corlocqt(_corlocqt, natm,8,3);
    Mat<int> loc_qt(_loc_qt, natm,8,3);
    for(int j=1; j<=charget._ny; ++j){
        charget(iatm,j) = charge[j-1];
        for (int k=1; k<=corlocqt._nz; ++k){
            corlocqt(iatm,j,k) = corlocq(j,k);
            loc_qt(iatm,j,k) = loc_q(j,k);

        }
    }

}

}//extern

    
