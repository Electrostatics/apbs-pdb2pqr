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
} Vcom;
extern Vcom* Vcom_ctor(int commtype);
extern int Vcom_size(Vcom *thee);
extern int Vcom_rank(Vcom *thee);

// Functions and Constructor from vmem.h:

typedef struct {
	Vmem();
} Vmem;
extern Vmem* Vmem_ctor(char *name);

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
extern int NOsh_parseFile(NOsh *thee, char *filename);

// Functions from routines.h:

typedef struct {
	AtomForce();
	~AtomForce();
} AtomForce;


// Functions for python implementation of objects that are arrays:
// Note: Currently does not support NOSH_MAXMOL, NOSH_MAXCALC
//	 Size specified in python file


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

void delete_Nosh(NOsh *nosh) {
    NOsh_dtor(&nosh);
}

void delete_Com(Vcom *com) {
    Vcom_dtor(&com);
}

void delete_Mem(Vmem *mem) {
    Vmem_dtor(&mem);
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

%inline %{

int parseInputFromString(NOsh *nosh, PyObject *string){

    int ret, bufsize;
    Vio *sock;
    
    startVio();
    bufsize = PyString_Size(string);

    VASSERT( bufsize <= VMAX_BUFSIZE );
    sock = Vio_ctor("BUFF","ASC",VNULL,"0","r");

    Vio_bufTake(sock, PyString_AsString(string), bufsize);

    ret = NOsh_parse(nosh, sock); 
    sock->VIObuffer = VNULL;
    Vio_dtor(&sock);
    return ret;
}

void Valist_load(Valist *thee, int size, PyObject *x, PyObject *y, PyObject *z, PyObject *chg, PyObject *rad){ 
   
    int i,j;
    double pos[3];

    Vatom *atom;
   
    VASSERT(thee != VNULL);

    thee->atoms = Vmem_malloc(thee->vmem, size, sizeof(Vatom));
    thee->number = size;
    for (i=0;i<size;i++){
        pos[0] = PyFloat_AsDouble(PyList_GetItem(x,i));
        pos[1] = PyFloat_AsDouble(PyList_GetItem(y,i));
        pos[2] = PyFloat_AsDouble(PyList_GetItem(z,i));
        Vatom_setCharge(&(thee->atoms[i]), PyFloat_AsDouble(PyList_GetItem(chg,i)));
        Vatom_setRadius(&(thee->atoms[i]), PyFloat_AsDouble(PyList_GetItem(rad,i)));
        Vatom_setPosition(&(thee->atoms[i]), pos);
        Vatom_setAtomID(&(thee->atoms[i]), i);
    }   

    thee->center[0] = 0.0;
    thee->center[1] = 0.0;
    thee->center[2] = 0.0;
    thee->maxrad = 0.0;
    thee->charge = 0.0;

    /* Reset stat variables */
    atom = &(thee->atoms[0]);
    for (i=0; i<3; i++) {
        thee->maxcrd[i] = thee->mincrd[i] = atom->position[i];
    }
    thee->maxrad = atom->radius;
   
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

int wrap_forceMG(Vmem *mem, NOsh *nosh, PBEparm *pbeparm, MGparm *mgparm,
 Vpmg *pmg, AtomForce *atomForce[NOSH_MAXCALC], Valist *alist[NOSH_MAXMOL], 
 int forcearray[NOSH_MAXCALC], int calcid)
{
    int *nforce;
    nforce = malloc(sizeof(int));
    *nforce = 0;
    forceMG(mem, nosh, pbeparm, mgparm, pmg, nforce, atomForce, alist);
    forcearray[calcid] = *nforce;

    return *nforce;
}

PyObject *getPotentials(NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg, Valist *alist){
    Vgrid *grid;
    Vatom *atom; 
    int i, rc, nx, ny, nz;
    double hx, hy, hzed, xcent, ycent, zcent, xmin, ymin, zmin;
    double value;
    double *position;
    PyObject *values;
    
    values = PyList_New(Valist_getNumberAtoms(alist));
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
        PyList_SetItem(values, i, PyFloat_FromDouble(value)); 
    } 
    Vgrid_dtor(&grid);    
    return values;
}

PyObject *getEnergies(Vpmg *pmg, Valist *alist){
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
//extern int forceMG(Vmem *mem, NOsh *nosh, PBEparm *pbeparm, MGparm *mgparm,
//  Vpmg *pmg, int *nforce, AtomForce *atomForce[NOSH_MAXCALC], Valist *alist[NOSH_MAXMOL]);
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
   int *INPUT, double *INOUT, double *INPUT, double *INPUT,
   double *INPUT);

%typemap(in) double [ANY] {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (double *) malloc((size+1)*sizeof(double));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyFloat_Check(o))
	    $1[i] = PyFloat_AsDouble(PyList_GetItem($input,i));
      else {
	    PyErr_SetString(PyExc_TypeError,"list must contain floats");
	    free($1);
	    return NULL;
      }
    }
    $1[i] = 0;
  } else {
        PyErr_SetString(PyExc_TypeError,"not a list");
        return NULL;
  }
}

extern int printEnergy(Vcom *com, NOsh *nosh, double totEnergy[NOSH_MAXCALC],
    int i);

