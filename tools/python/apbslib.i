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

// Functions from routines.h:

typedef struct {
	AtomForce();
	~AtomForce();
} AtomForce;


// Functions for python implementation of objects that are arrays:
// Note: Currently does not support NOSH_MAXMOL, NOSH_MAXCALC
//	 Size specified in python file

%include pointer.i

%inline %{
Valist **new_valist(int maxargs) {
   return (Valist **) malloc(maxargs*sizeof(Valist *));
}

Vgrid **new_gridlist(int maxargs) {
   return (Vgrid **) malloc(maxargs*sizeof(Vgrid *));
}

Vpmg **new_pmglist(int maxargs) {
   return (Vpmg **) malloc(maxargs*sizeof(Vpmg *));
}

Vpmg *get_Vpmg(Vpmg **args, int n) {
   return (Vpmg *)args[n];
}

Vpmgp **new_pmgplist(int maxargs) {
   return (Vpmgp **) malloc(maxargs*sizeof(Vpmgp *));
}

Vpbe **new_pbelist(int maxargs) {
   return (Vpbe **) malloc(maxargs*sizeof(Vpbe *));
}

Vpbe *get_Vpbe(Vpbe **args, int n) { 
    return (Vpbe *)args[n];
}

AtomForce **new_atomforcelist(int maxargs) {
   return (AtomForce **) malloc(maxargs*sizeof(AtomForce *));
}

void delete_atomforcelist(AtomForce **a) {
    free(a);
  } 
void delete_valist(Valist **a) {
    free(a);
  } 
void delete_gridlist(Vgrid **a) {
    free(a);
  } 
void delete_pmglist(Vpmg **a) {
    free(a);
  } 
void delete_pmgplist(Vpmgp **a) {
    free(a);
  }
void delete_pbelist(Vpbe **a) {
    free(a);
  }

AtomForce **get_AtomForce(AtomForce **aforce, int n){
    return &aforce[n];
}

Valist *make_Valist(Valist **args, int n){
    args[n] = Valist_ctor();    
    return args[n];
}

// Generic array of doubles and ints:
//   Constructors, Destructors, Gets, and Sets


double *double_array(int size) {
     return (double *) malloc(size*sizeof(double));
  }

int *int_array(int size){
     return (int *) malloc(size*sizeof(int));
  }

void delete_double_array(double *d) {
    free(d);
  } 

void delete_int_array(int *i) {
    free(i);
  } 

double get_entry(double *array, int i){
	    return array[i];
  }

void set_entry(double *array, int i, double val){
	    array[i] = val;
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
extern int npenergyMG(NOsh* nosh, int icalc, Vpmg *pmg, int *nenergy, double *npEnergy);
extern void killEnergy();
extern int forceMG(Vmem *mem, NOsh *nosh, PBEparm *pbeparm, MGparm *mgparm,
  Vpmg *pmg, int *nforce, AtomForce *atomForce[NOSH_MAXCALC], Valist *alist[NOSH_MAXMOL]);
extern void killForce(Vmem *mem, NOsh *nosh, int nforce[NOSH_MAXCALC],
  AtomForce *atomForce[NOSH_MAXCALC]);
extern int writedataMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg);
extern int writematMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg);
extern int printForce(Vcom *com, NOsh *nosh, int nforce[NOSH_MAXCALC],
  AtomForce *atomForce[NOSH_MAXCALC], int i);
extern void startVio();
extern double Vacc_molAcc(Vacc *thee, double center[3], double radius);
extern double Vacc_vdwAcc(Vacc *thee, double center[3]);


// Typemaps and functions for easy Python Access

%include typemaps.i

extern int energyMG(NOsh* nosh, int icalc, Vpmg *pmg,
   int *INPUT, double *BOTH, double *INPUT, double *INPUT,
   double *INPUT);

%typemap(python,in) double * {
  /* Check if is a list */
  if (PyList_Check($source)) {
    int size = PyList_Size($source);
    int i = 0;
    $target = (double *) malloc((size+1)*sizeof(double));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($source,i);
      if (PyFloat_Check(o))
	    $target[i] = PyFloat_AsDouble(PyList_GetItem($source,i));
      else {
	    PyErr_SetString(PyExc_TypeError,"list must contain floats");
	    free($target);
	    return NULL;
      }
    }
    $target[i] = 0;
  } else {
        PyErr_SetString(PyExc_TypeError,"not a list");
        return NULL;
  }
}

extern int printEnergy(Vcom *com, NOsh *nosh, double totEnergy[NOSH_MAXCALC],
  int i);

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

PyObject *getPotentials(NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg, Valist *alist){
    Vatom *atom; 
    int i;
    double energy;
    PyObject *values;
    
    values = PyList_New(Valist_getNumberAtoms(alist)); 
    for (i=0;i<Valist_getNumberAtoms(alist);i++){
        atom = Valist_getAtom(alist, i);
        energy = Vpmg_qfAtomEnergy(pmg, atom);
        PyList_SetItem(values, i, PyFloat_FromDouble(energy));
    } 
    return values;
}

PyObject *getForces(AtomForce **atomForce, Valist *alist){
    int i, j;
    PyObject *dict;
    PyObject *qfvalues, *qf, *qfholder;
    PyObject *ibvalues, *ib, *ibholder;
    PyObject *dbvalues, *db, *dbholder;
    PyObject *npvalues, *np, *npholder;
    
    dict = PyDict_New(); 
    qfvalues = PyList_New(Valist_getNumberAtoms(alist));
    dbvalues = PyList_New(Valist_getNumberAtoms(alist)); 
    ibvalues = PyList_New(Valist_getNumberAtoms(alist)); 
    npvalues = PyList_New(Valist_getNumberAtoms(alist)); 
    
    qfholder = PyList_New(3);
    dbholder = PyList_New(3);
    ibholder = PyList_New(3);
    npholder = PyList_New(3);
    
    qf = PyString_FromString("qf");
    db = PyString_FromString("db");
    ib = PyString_FromString("ib");
    np = PyString_FromString("np");

    for (i=0;i<Valist_getNumberAtoms(alist);i++){
        for (j=0;j<3;j++){
           PyList_SetItem(qfholder, j, PyFloat_FromDouble(atomForce[0][i].qfForce[j]));
           PyList_SetItem(dbholder, j, PyFloat_FromDouble(atomForce[0][i].dbForce[j]));
           PyList_SetItem(ibholder, j, PyFloat_FromDouble(atomForce[0][i].ibForce[j]));
           PyList_SetItem(npholder, j, PyFloat_FromDouble(atomForce[0][i].npForce[j]));
        }
        PyList_SetItem(qfvalues, i,  PyList_GetSlice(qfholder, 0, 3));
        PyList_SetItem(dbvalues, i,  PyList_GetSlice(dbholder, 0, 3));
        PyList_SetItem(ibvalues, i,  PyList_GetSlice(ibholder, 0, 3));
        PyList_SetItem(npvalues, i,  PyList_GetSlice(npholder, 0, 3));
    }
    PyDict_SetItem(dict, qf, qfvalues);
    PyDict_SetItem(dict, db, dbvalues);
    PyDict_SetItem(dict, np, npvalues);
    PyDict_SetItem(dict, ib, ibvalues);
    return dict;
}
%}
