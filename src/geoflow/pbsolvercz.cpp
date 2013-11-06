#include "modules.h"
#include "Mat.h"
#include <valarray>
#include <vector>

#include <iostream>
#include <fstream>
#include <Eigen/Sparse>


extern "C"{

void pbsolver(double* _eps, double* _phi, double* bgf, int nx,int ny,int nz, double dcel, double tol, int iter){
    int n = nx*ny*nz;
    std::vector<double> _eps1(n), _eps2(n), _eps3(n);
    Mat<> eps(_eps,nx,ny,nz), eps1(_eps1,nx,ny,nz), eps2(_eps2,nx,ny,nz), eps3(_eps3,nx,ny,nz);
    for(int i=1; i<nx; ++i){
    for(int j=1; j<ny; ++j){
    for(int k=1; k<nz; ++k){
        eps1(i,j,k) = (eps(i+1,j,k) + eps(i,j,k))/2.0;
        eps2(i,j,k) = (eps(i,j+1,k) + eps(i,j,k))/2.0;
        eps3(i,j,k) = (eps(i,j,k+1) + eps(i,j,k))/2.0;
    }}}

    std::vector< Eigen::Triplet<double> > tripletList;
    tripletList.reserve(nx*ny*nz);
    Eigen::VectorXd b(nx*ny*nz);
    Mat<> phi(_phi, nx,ny,nz);
    Eigen::VectorXd phi_flat(nx*ny*nz);

    for(int i=0; i<nx; ++i){
    for(int j=0; j<ny; ++j){
    for(int k=0; k<nz; ++k){
        int ijk = i*nz*ny + j*nz + k;
        if(i==0 || i==nx-1 || j==0 || j==ny-1 || k==0 || k==nz-1){
            tripletList.push_back( Eigen::Triplet<double>(ijk, ijk, 1.0) );
        }else{
            double foo = -(  eps1(i+1,j+1,k+1) + eps1(i,j+1,k+1)
                           + eps2(i+1,j+1,k+1) + eps2(i+1,j,k+1)
                           + eps3(i+1,j+1,k+1) + eps3(i+1,j+1,k) )/dcel/dcel;
            tripletList.push_back( Eigen::Triplet<double>(ijk, ijk, foo) ); 

            std::valarray<double> weit(6);
            weit[0] = eps1(i,j+1,k+1);
            weit[1] = eps2(i+1,j,k+1);
            weit[2] = eps3(i+1,j+1,k);
            weit[3] = eps3(i+1,j+1,k+1);
            weit[4] = eps2(i+1,j+1,k+1);
            weit[5] = eps1(i+1,j+1,k+1);
            weit /= dcel*dcel;

            int jj = ijk - nz*ny; 
            if(jj>=0){ tripletList.push_back( Eigen::Triplet<double>(ijk, jj, weit[0]) ); } 
            
            jj = ijk - nz; 
            if(jj>=0){ tripletList.push_back( Eigen::Triplet<double>(ijk, jj, weit[1]) ); } 

            jj = ijk - 1;
            if(jj>=0){ tripletList.push_back( Eigen::Triplet<double>(ijk, jj, weit[2]) ); } 

            jj = ijk + 1; 
            if(jj<n){ tripletList.push_back( Eigen::Triplet<double>(ijk, jj, weit[3]) ); } 

            jj = ijk + nz; 
            if(jj<n){ tripletList.push_back( Eigen::Triplet<double>(ijk, jj, weit[4]) ); } 
            
            jj = ijk + nz*ny; 
            if(jj<n){ tripletList.push_back( Eigen::Triplet<double>(ijk, jj, weit[5]) ); } 

        }
        b(ijk) = bgf[ijk];
        phi_flat(ijk) = phi(i+1,j+1,k+1);
    }}}

    Eigen::SparseMatrix<double> A(n, n);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    A.makeCompressed();

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver(A);
    solver.setMaxIterations(iter);
    solver.setTolerance(tol);

    phi_flat = solver.solveWithGuess(b, phi_flat);

    for(int i=1; i<=nx; ++i){
    for(int j=1; j<=ny; ++j){
    for(int k=1; k<=nz; ++k){
        int ijk = (i-1)*nz*ny + (j-1)*nz + k;
        phi(i,j,k) = phi_flat(ijk-1);
    }}}
}

}//extern
