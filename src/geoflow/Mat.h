///  @file    Mat.h
///  @author  Andrew Stevens
///  @brief Wrapper class for 2d and 3d arrays
///  @ingroup Geoflow
///  @version $Id$
///  @attention
///  @verbatim
///
/// APBS -- Adaptive Poisson-Boltzmann Solver
///
///  Nathan A. Baker (nathan.baker@pnnl.gov)
///  Pacific Northwest National Laboratory
///
///  Additional contributing authors listed in the code documentation.
///
/// Copyright (c) 2010-2012 Battelle Memorial Institute. Developed at the
/// Pacific Northwest National Laboratory, operated by Battelle Memorial
/// Institute, Pacific Northwest Division for the U.S. Department of Energy.
///
/// Portions Copyright (c) 2002-2010, Washington University in St. Louis.
/// Portions Copyright (c) 2002-2010, Nathan A. Baker.
/// Portions Copyright (c) 1999-2002, The Regents of the University of
/// California.
/// Portions Copyright (c) 1995, Michael Holst.
/// All rights reserved.
///
/// Redistribution and use in source and binary forms, with or without
/// modification, are permitted provided that the following conditions are met:
///
/// Redistributions of source code must retain the above copyright notice, this
/// list of conditions and the following disclaimer.
///
/// Redistributions in binary form must reproduce the above copyright notice,
/// this list of conditions and the following disclaimer in the documentation
/// and/or other materials provided with the distribution.
///
/// Neither the name of the developer nor the names of its contributors may be
/// used to endorse or promote products derived from this software without
/// specific prior written permission.
///
/// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
/// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
/// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
/// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
/// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
/// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
/// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
/// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
/// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
/// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
/// THE POSSIBILITY OF SUCH DAMAGE.
///
/// @endverbatim

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

