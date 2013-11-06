#ifndef MAT_H
#define MAT_H

#include <valarray>
#include <vector>
#include <cmath>

template<typename T=double> struct Mat;
template<typename T> struct Mat{
    int _nx,_ny,_nz;
    T* _data;
    Mat(T* data, int nx, int ny, int nz=0):  _nx(nx), _ny(ny), _nz(nz), _data(data) {}
    Mat(std::vector<T>& val, int nx, int ny, int nz=0): _nx(nx), _ny(ny), _nz(nz), _data(val.data()) {}

    int index(int x, int y, int z) const { return x-1 + (y-1)*_nx + (z-1)*_nx*_ny; }
    T& operator()(int x, int y, int z=1){ return _data[index(x,y,z)]; }
    T operator()(int x, int y, int z=1) const { return _data[index(x,y,z)]; }
    
    T& operator[](int i){ return _data[i]; }
    T operator[](int i) const { return _data[i]; }

    T* data(){return _data;}
    int size(){return _nx*_ny*(_nz==0?1:_nz);}
    T* end(){return _data + size();}

    void operator=(T a){ std::fill(data(), end(), a); }
    void operator=(Mat& a){ std::copy(a.data(), a.end(), data()); }
};

#endif

