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
} NOsh;

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
