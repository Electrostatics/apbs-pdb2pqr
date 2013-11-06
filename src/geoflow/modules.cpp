#include "modules.h"

Comdata comdata;
LJ lj;

double dot(double x, double y, double z){
    return x*x + y*y + z*z;
}

double xvalue(int& i){ return (i - 1)*comdata.deltax + comdata.xleft; }
double yvalue(int& i){ return (i - 1)*comdata.deltay + comdata.yleft; }
double zvalue(int& i){ return (i - 1)*comdata.deltaz + comdata.zleft; }

