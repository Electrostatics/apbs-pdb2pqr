/* Input file for creating Python wrappers for APBS via swig
   Author: Todd Dolinsky
   Email: todd@ccb.wustl.edu
*/

/* 
Header files:
-----------------------
*/ 
     
%{
#define APBS_SWIG 1
#include "maloc/maloc.h"
#include "apbscfg.h" 
#include "routines.h"
%} 

/* 
Structures and functions to be wrapped:
--------------------------------------
*/
   
#define APBS_SWIG 1
#define VEXTERNC extern  

// Functions and Constructors from mgparm.h:
  
typedef struct {
	MGparm();
	~MGparm();
	MGparm_CalcType type;                                           
} MGparm; 
      
// Functions and Constructors from pbeparm.h:
  
typedef struct {  
	PBEparm();
	~PBEparm();
	double temp;
} PBEparm;

// Functions and Constructor from vcom.h:

typedef struct {
	Vcom();
	~Vcom();
} Vcom;
extern Vcom* Vcom_ctor(int commtype);
extern void Vcom_dtor(Vcom **thee);
extern int Vcom_size(Vcom *thee);
extern int Vcom_rank(Vcom *thee);

// Functions and Constructor from vmem.h:

typedef struct {
	Vmem();
	~Vmem();
} Vmem;
extern Vmem* Vmem_ctor(char *name);
extern void Vmem_dtor(Vmem **thee);

// Functions and Constructor from vpmg.h:

typedef struct {
	Vpmg();
	~Vpmg();
} Vpmg;

// Functions and Constructor from vpbe.h:
typedef struct {
	Vpbe();
	~Vpbe();
    Vacc *acc;
} Vpbe;

// Functions and Constructors from nosh.h:

typedef struct { 
	MGparm *mgparm;         
	FEMparm *femparm;       
	PBEparm *pbeparm;       
	int calctype;           
} NOsh_calc;

typedef struct {
	NOsh();
	~NOsh();
	int ncalc;
	int nprint;             
    int nelec;
    NOsh_PrintType printwhat[NOSH_MAXPRINT];
} NOsh;

enum NOsh_PrintType {
    NPT_ENERGY=0, NPT_FORCE=1
};


#if !defined(VINLINE_NOSH)
extern NOsh_calc* NOsh_getCalc(NOsh *thee, int icalc);
#else
#define NOsh_getCalc(thee, icalc) ((thee)->calc[(icalc)])
#endif

extern char* NOsh_elecname(NOsh *thee, int ielec);
extern int NOsh_elec2calc(NOsh *thee, int icalc);
extern NOsh_PrintType NOsh_printWhat(NOsh *thee, int iprint); 
extern int NOsh_ctor2(NOsh *thee, int rank, int size);
extern void NOsh_dtor(NOsh **thee); 
extern int NOsh_parseFile(NOsh *thee, char *filename);

// Functions for python implementation of objects that are arrays:
// Note: Currently does not support NOSH_MAXMOL, NOSH_MAXCALC
//	 Size specified in python file

%include pointer.i

%inline %{
Valist **new_valist(int maxargs) {
   return (Valist **) malloc(maxargs*sizeof(Valist *));
}
%}

%inline %{
Vgrid **new_gridlist(int maxargs) {
   return (Vgrid **) malloc(maxargs*sizeof(Vgrid *));
}
%}


%inline %{
Vpmg **new_pmglist(int maxargs) {
   return (Vpmg **) malloc(maxargs*sizeof(Vpmg *));
}

Vpmg *get_Vpmg(Vpmg **args, int n) {
   return (Vpmg *)args[n];
}
%}

%inline %{
Vpmgp **new_pmgplist(int maxargs) {
   return (Vpmgp **) malloc(maxargs*sizeof(Vpmgp *));
}
%}

%inline %{
Vpbe **new_pbelist(int maxargs) {
   return (Vpbe **) malloc(maxargs*sizeof(Vpbe *));
}
%}

%inline %{
Vpbe *get_Vpbe(Vpbe **args, int n) { 
    return (Vpbe *)args[n];
}
%}



%inline %{
AtomForce **new_atomforcelist(int maxargs) {
   return (AtomForce **) malloc(maxargs*sizeof(AtomForce *));
}
%}

// Stubs for deleting structures - not currently used
/*
%inline %{
NOsh **getNosh(NOsh *thee){
	return &thee;
}
%}
/*
%inline %{
Vcom **getCom(Vcom *thee){
   return &thee;
}
%}

%inline %{
Vmem **getMem(Vmem *thee){
   return &thee;
}
%}
*/

// Generic array of doubles:

%inline %{
double *double_array(int size) {
     return (double *) malloc(size*sizeof(double));
  }
%}

// Generic array of ints:

%inline %{
int *int_array(int size){
     return (int *) malloc(size*sizeof(int));
}
%}

// Functions for PDB2PQR interface

%inline %{
void Valist_load(Valist *thee, int size, double *x, double *y, double *z, double *chg, double *rad){ 
    
    Vatom *atoms = VNULL;
    Vatom *atom;
    int i;
    int j;

    double pos[3];
    atoms = Vmem_malloc(thee->vmem, size, sizeof(Vatom));
    thee->number = 0;
    for (i=0;i<size;i++){
        pos[0] = x[i];
        pos[1] = y[i];
        pos[2] = z[i];
        Vatom_setCharge(&(atoms[thee->number]), chg[i]);
        Vatom_setRadius(&(atoms[thee->number]), rad[i]);
        Vatom_setPosition(&(atoms[thee->number]), pos);
        (thee->number)++;
      
    }
  
    thee->atoms = Vmem_malloc(thee->vmem, thee->number,(sizeof(Vatom)));
    VASSERT(thee->atoms != VNULL);
    for (i=0; i<thee->number; i++) {
        Vatom_copyTo(&(atoms[i]), &(thee->atoms[i]));
        Vatom_dtor2(&(atoms[i]));
    }
    Vmem_free(thee->vmem, size, sizeof(Vatom), (void **)&atoms);
    
    VASSERT(thee != VNULL);

    thee->center[0] = 0.;
    thee->center[1] = 0.;
    thee->center[2] = 0.;
    thee->maxrad = 0.;
    thee->charge = 0.;

    /* Reset stat variables */
    atom = &(thee->atoms[0]);
    for (i=0; i<3; i++) {
        thee->maxcrd[i] = thee->mincrd[i] = atom->position[i];
    }
    thee->maxrad = atom->radius;
    thee->charge = 0.0;
   
    for (i=0; i<thee->number; i++) {

        atom = &(thee->atoms[i]);
        for (j=0; j<3; j++) {
            if (atom->position[j] < thee->mincrd[j]) 
              thee->mincrd[j] = atom->position[j];
            if (atom->position[j] > thee->maxcrd[j]) 
              thee->maxcrd[j] = atom->position[j];
        }
        if (atom->radius > thee->maxrad) thee->maxrad = atom->radius;
        thee->charge = thee->charge + atom->charge;
    } 
  
    thee->center[0] = 0.5*(thee->maxcrd[0] + thee->mincrd[0]);
    thee->center[1] = 0.5*(thee->maxcrd[1] + thee->mincrd[1]);
    thee->center[2] = 0.5*(thee->maxcrd[2] + thee->mincrd[2]);
}
%}

%inline %{
double *getPotentials(NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg, Valist *alist){
    Vgrid *grid;
    Vatom *atom; 
    int i, rc, nx, ny, nz;
    double hx, hy, hzed, xcent, ycent, zcent, xmin, ymin, zmin;
    double value;
    double *position, *values;
    
    values = Vmem_malloc(alist->vmem, Valist_getNumberAtoms(alist),(sizeof(double)));
    nx = pmg->pmgp->nx;
    ny = pmg->pmgp->ny;
    nz = pmg->pmgp->nz;
    hx = pmg->pmgp->hx;
    hy = pmg->pmgp->hy;
    hzed = pmg->pmgp->hzed;
    xcent = pmg->pmgp->xcent;
    ycent = pmg->pmgp->ycent;
    zcent = pmg->pmgp->zcent;
    xmin = xcent - 0.5*(nx-1)*hx;
    ymin = ycent - 0.5*(ny-1)*hy;
    zmin = zcent - 0.5*(nz-1)*hzed;
   
    Vpmg_fillArray(pmg, pmg->rwork, VDT_POT, 0.0, pbeparm->pbetype);
    grid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin,
                  pmg->rwork);
    for (i=0;i<Valist_getNumberAtoms(alist);i++){
        atom = Valist_getAtom(alist, i);
        position = Vatom_getPosition(atom); 
        Vgrid_value(grid, position, &value);
        values[i] = value; 
    } 
    Vgrid_dtor(&grid);    
    return values;
}
%}

%inline %{
double **getqfForces(AtomForce **atomForce, Valist *alist){
    int i, j;
    double **values;
    double *holder; 
  
    holder = Vmem_malloc(alist->vmem, 3, sizeof(double));
    values = Vmem_malloc(alist->vmem, Valist_getNumberAtoms(alist),(sizeof(holder)));
    for (i=0;i<Valist_getNumberAtoms(alist);i++){
        for (j=0;j<3;j++){
            holder[j] = (*atomForce)[i].qfForce[j];
        }
        values[i] = holder;
    }   
    return values;
}
%}

%inline %{
double **getibForces(AtomForce **atomForce, Valist *alist){
    int i, j;
    double **values;
    double *holder; 
  
    holder = Vmem_malloc(alist->vmem, 3, sizeof(double));
    values = Vmem_malloc(alist->vmem, Valist_getNumberAtoms(alist),(sizeof(holder)));
    for (i=0;i<Valist_getNumberAtoms(alist);i++){
        for (j=0;j<3;j++){
            holder[j] = (*atomForce)[i].ibForce[j];
        }
        values[i] = holder;
    }   
    return values;
}
%}

%inline %{
double **getdbForces(AtomForce **atomForce, Valist *alist){
    int i, j;
    double **values;
    double *holder; 
  
    holder = Vmem_malloc(alist->vmem, 3, sizeof(double));
    values = Vmem_malloc(alist->vmem, Valist_getNumberAtoms(alist),(sizeof(holder)));
    for (i=0;i<Valist_getNumberAtoms(alist);i++){
        for (j=0;j<3;j++){
            holder[j] = (*atomForce)[i].dbForce[j];
        }
        values[i] = holder;
    }   
    return values;
}
%}

%inline %{
double **getnpForces(AtomForce **atomForce, Valist *alist){
    int i, j;
    double **values;
    double *holder; 
  
    holder = Vmem_malloc(alist->vmem, 3, sizeof(double));
    values = Vmem_malloc(alist->vmem, Valist_getNumberAtoms(alist),(sizeof(holder)));
    for (i=0;i<Valist_getNumberAtoms(alist);i++){
        for (j=0;j<3;j++){
            holder[j] = (*atomForce)[i].npForce[j];
        }
        values[i] = holder;
    }   
    return values;
}
%}


%inline %{
double *get_double_entry(double **array, int i){
	    return array[i];
  }
%}

%inline %{
double get_entry(double *array, int i){
	    return array[i];
  }
%}

%inline %{
void set_entry(double *array, int i, double val){
	    array[i] = val;
  }
%}

%inline %{
Valist *make_Valist(Valist **args, int n){
    args[n] = Valist_ctor();    
    return args[n];
}
%}

// Additional functions for reading input from buffers

extern int NOsh_parse(NOsh *thee, Vio *sock);

%inline %{
Vio * Vio_setup(char *key, const char *iodev, const char *iofmt, const char *iohost, const char *iofile, char * string){
    Vio *sock = VNULL;
    char buf[VMAX_BUFSIZE];
    int bufsize = 0;
    bufsize = strlen(string);
    VASSERT( bufsize <= VMAX_BUFSIZE );
    strncpy(buf, string, VMAX_BUFSIZE);
    VASSERT( VNULL != (sock=Vio_socketOpen(key,iodev,iofmt,iohost,iofile)));
    Vio_bufTake(sock, buf, bufsize);
    return sock;
}
%}

// Functions from routines.h:

typedef struct {
	AtomForce();
	~AtomForce();
} AtomForce;

extern int loadMolecules(NOsh *nosh, Valist *alist[NOSH_MAXMOL]);
extern void killMolecules(NOsh *nosh, Valist *alist[NOSH_MAXMOL]);
extern int loadDielMaps(NOsh *nosh, Vgrid *dielXMap[NOSH_MAXMOL],
Vgrid *dielYMap[NOSH_MAXMOL], Vgrid *dielZMap[NOSH_MAXMOL]);
extern void killDielMaps(NOsh *nosh, Vgrid *dielXMap[NOSH_MAXMOL],
Vgrid *dielYMap[NOSH_MAXMOL], Vgrid *dielZMap[NOSH_MAXMOL]);
extern int loadKappaMaps(NOsh *nosh, Vgrid *kappa[NOSH_MAXMOL]);
extern void killKappaMaps(NOsh *nosh, Vgrid *kappa[NOSH_MAXMOL]);
extern int loadChargeMaps(NOsh *nosh, Vgrid *charge[NOSH_MAXMOL]);
extern void killChargeMaps(NOsh *nosh, Vgrid *charge[NOSH_MAXMOL]);
extern void printPBEPARM(PBEparm *pbeparm);
extern void printMGPARM(MGparm *mgparm, double realCenter[3]);
extern int initMG(int i, NOsh *nosh, MGparm *mgparm,
  PBEparm *pbeparm, double realCenter[3], Vpbe *pbe[NOSH_MAXCALC],
  Valist *alist[NOSH_MAXMOL], Vgrid *dielXMap[NOSH_MAXMOL], 
  Vgrid *dielYMap[NOSH_MAXMOL], Vgrid *dielZMap[NOSH_MAXMOL], 
  Vgrid *kappaMap[NOSH_MAXMOL], Vgrid *chargeMap[NOSH_MAXMOL], 
  Vpmgp *pmgp[NOSH_MAXCALC], Vpmg *pmg[NOSH_MAXCALC]);
extern void killMG(NOsh *nosh, Vpbe *pbe[NOSH_MAXCALC],
  Vpmgp *pmgp[NOSH_MAXCALC], Vpmg *pmg[NOSH_MAXCALC]);
extern int solveMG(NOsh *nosh, Vpmg *pmg, MGparm_CalcType type);
extern int setPartMG(NOsh *nosh, MGparm *mgparm, Vpmg *pmg);
extern int energyMG(NOsh* nosh, int icalc, Vpmg *pmg,
  int *nenergy, double *totEnergy, double *qfEnergy, double *qmEnergy,
  double *dielEnergy);
extern int npenergyMG(NOsh* nosh, int icalc, Vpmg *pmg, int *nenergy, double *npEnergy);
extern void killEnergy();
extern int forceMG(Vmem *mem, NOsh *nosh, PBEparm *pbeparm, MGparm *mgparm,
  Vpmg *pmg, int *nforce, AtomForce **atomForce, Valist *alist[NOSH_MAXMOL]);
extern void killForce(Vmem *mem, NOsh *nosh, int nforce[NOSH_MAXCALC],
  AtomForce *atomForce[NOSH_MAXCALC]);
extern int writedataMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg);
extern int writematMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg);
extern int printEnergy(Vcom *com, NOsh *nosh, double totEnergy[NOSH_MAXCALC],
  int i);
extern int printForce(Vcom *com, NOsh *nosh, int nforce[NOSH_MAXCALC],
  AtomForce *atomForce[NOSH_MAXCALC], int i);
extern void startVio();
extern double Vacc_molAcc(Vacc *thee, double center[3], double radius);
extern double Vacc_vdwAcc(Vacc *thee, double center[3]);
