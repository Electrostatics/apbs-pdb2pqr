#ifndef MODULES_H
#define MODULES_H

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
    double* xc, *yc, *zc;
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
void domainini(double xyzr[MAXATOMS][XYZRWIDTH], const int natm, const double extvalue);
void chargedist(double* atmpos, double* chratm, int& natm, double* charget, double* corlocqt, int* loc_qt, int& iatm);
void yhsurface(double xyzr[MAXATOMS][XYZRWIDTH], double* ljepsilon, int natm, double tott,
    double deltat, double* phix, double* surfu, int i, double& area, double& vol, double& attint,
    double alpha, int iadi, int igfin);
void seteqb(double* bg, double xyzr[MAXATOMS][XYZRWIDTH], double* pqr, int& natm, double* charget, double* corlocqt, double *epsilonsp);
//  void writerms_gama( double* sumpot, double* expv, double* elec, int* natm, double* gama, int *ngiter, double *pres );
void pbsolver( double* eps, double* phi, double* bg, int nx, int ny, int nz, double dcel,  double tol, int iter);

double xvalue(int& i);
double yvalue(int& i);
double zvalue(int& i);

int inverx(double& x);
int invery(double& y);
int inverz(double& z);
}

#endif

