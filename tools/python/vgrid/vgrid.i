/* Input file for creating Python wrappers for Vgrid via swig
   Author: Todd Dolinsky
   Email: todd@ccb.wustl.edu
 
Header files:
-----------------------
*/ 

%{
#include "routines.h"
#include "apbs/vgrid.h"
%}

#define VEXTERNC extern

%include pointer.i

// Generic array of doubles:

%inline %{
double *null_array(){
     return NULL;
}
%}

%inline %{
double *double_array(int size) {
     return (double *) malloc(size*sizeof(double));
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

typedef struct {
	Vgrid();
	~Vgrid();
	int nx;
	int ny;
	int nz;
	double hx;
	double hy;
	double hzed;
	double xmin;
	double ymin;
	double zmin;
    double *data;
} Vgrid;

extern Vgrid* Vgrid_ctor(int nx, int ny, int nz, double hx, double hy,
						 double hzed, double xmin, double ymin, double zmin, 
						 double *data);
extern int Vgrid_ctor2(Vgrid *thee, int nx, int ny, int nz, double hx, 
					double hy, double hzed, double xmin, double ymin, 
					double zmin, double *data);
extern int Vgrid_value(Vgrid *thee, double x[3], double *value);
extern void Vgrid_dtor(Vgrid **thee);
extern void Vgrid_dtor2(Vgrid *thee);
extern int Vgrid_curvature(Vgrid *thee, double pt[3], int cflag, double *curv);
extern int Vgrid_gradient(Vgrid *thee, double pt[3], double grad[3]);

extern void Vgrid_writeUHBD(Vgrid *thee, const char *iodev, const char *iofmt, const char *thost, const char *fname, char *title, double *pvec);

extern void Vgrid_writeDX(Vgrid *thee, const char *iodev, const char *iofmt, const char *thost, const char *fname, char *title, double *pvec);

extern int Vgrid_readDX(Vgrid *thee, const char *iodev, const char *iofmt, const char *thost, const char *fname);

extern void startVio();
