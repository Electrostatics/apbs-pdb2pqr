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
    size_t _nx,_ny,_nz;
    T* _data;
    Mat(T* data, size_t nx, size_t ny, size_t nz=1):  _nx(nx), _ny(ny), _nz(nz), _data(data) {}
    Mat(std::vector<T>& val, size_t nx, size_t ny, size_t nz=0): _nx(nx), _ny(ny), _nz(nz), _data(val.data()) {}

    size_t index(size_t x, size_t y, size_t z) const { return x-1 + _nx*((y-1) + _ny*(z-1)); }
    T& operator()(size_t x, size_t y, size_t z=1){ return _data[index(x,y,z)]; }
    T operator()(size_t x, size_t y, size_t z=1) const { return _data[index(x,y,z)]; }
    
    T& operator[](size_t i){ return _data[i]; }
    T operator[](size_t i) const { return _data[i]; }

    T* data(){return _data;}
    size_t size(){return _nx*_ny*_nz;}
    T* end(){return _data + size();}

    void operator=(T a){ std::fill(data(), end(), a); }
    void operator=(Mat& a){ std::copy(a.data(), a.end(), data()); }
};

template<typename T> struct StencilBase{
    T *c, *xm,*ym,*zm, *xp,*yp,*zp,
      *xym,*xzm,*yzm, *xyp,*xzp,*yzp,
      *xm_yp,*xm_zp,*ym_zp, *xp_ym,*xp_zm,*yp_zm;
    
}

template<typename T> struct Stencil:  public std::iterator<std::forward_iterator_tag, Mat<T>>
{
    Mat<T>& _mat;
    T h,halfh,h2,qrth2;
    size_t i;

    T *c, *xm,*ym,*zm, *xp,*yp,*zp,
      *xym,*xzm,*yzm, *xyp,*xzp,*yzp,
      *xm_yp,*xm_zp,*ym_zp, *xp_ym,*xp_zm,*yp_zm;

    Stencil(Mat<T>& mat, T _h): _mat(mat), h(1.0/_h), halfh(1.0/(2*_h)), h2(1./(_h*_h)),
        qrth2(1.0/(4*_h*_h)), i(_mat.index(2,2,2)),
        c(_mat.data() + _mat.index(2,2,2)),
        xm(_mat.data() + _mat.index(1,2,2)), ym(_mat.data() + _mat.index(2,1,2)), zm(_mat.data() + _mat.index(2,2,1)),
        xp(_mat.data() + _mat.index(3,2,2)), yp(_mat.data() + _mat.index(2,3,2)), zp(_mat.data() + _mat.index(2,2,3)),
        xym(_mat.data() + _mat.index(1,1,2)), xzm(_mat.data() + _mat.index(1,2,1)), yzm(_mat.data() + _mat.index(2,1,1)),
        xyp(_mat.data() + _mat.index(3,3,2)), xzp(_mat.data() + _mat.index(3,2,3)), yzp(_mat.data() + _mat.index(2,3,3)),
        xm_yp(_mat.data() + _mat.index(1,3,2)), xm_zp(_mat.data() + _mat.index(1,2,3)), ym_zp(_mat.data() + _mat.index(2,1,3)),
        xp_ym(_mat.data() + _mat.index(3,1,2)), xp_zm(_mat.data() + _mat.index(3,2,1)), yp_zm(_mat.data() + _mat.index(2,3,1))
    {}

    void next(){
        ptrStep(1);
        if((i+1)%_mat._nx == 0){
            if( (i+1+_mat._nx)%(_mat._nx*_mat._ny) == 0 ){
                ptrStep(_mat._nx + 2);
            }else{
                ptrStep(2);
            }
        }
        
//         c = _mat(x,y,z),
//         xm = _mat(x-1,y,z); ym = _mat(x,y-1,z); zm = _mat(x,y,z-1);
//         xp = _mat(x+1,y,z); yp = _mat(x,y+1,z); zp = _mat(x,y,z+1);
//         xym = _mat(x-1,y-1,z); xzm = _mat(x-1,y,z-1); yzm = _mat(x,y-1,z-1);
//         xyp = _mat(x+1,y+1,z); xzp = _mat(x+1,y,z+1); yzp = _mat(x,y+1,z+1);
//         xm_yp = _mat(x-1,y+1,z); xm_zp = _mat(x-1,y,z+1); ym_zp = _mat(x,y-1,z+1);
//         xp_ym = _mat(x+1,y-1,z); xp_zm = _mat(x+1,y,z-1); yp_zm = _mat(x,y+1,z-1);
    }

    void ptrStep(size_t numSteps){
        i += numSteps;
        c += numSteps;
        xm += numSteps; ym += numSteps; zm += numSteps;
        xp += numSteps; yp += numSteps; zp += numSteps;
        xym += numSteps; xzm += numSteps; yzm += numSteps;
        xyp += numSteps; xzp += numSteps; yzp += numSteps;
        xm_yp += numSteps; xm_zp += numSteps; ym_zp += numSteps;
        xp_ym += numSteps; xp_zm += numSteps; yp_zm += numSteps;
    }

    T dx() { return halfh*(*xp - *xm); }
    T dy() { return halfh*(*yp - *ym); }
    T dz() { return halfh*(*zp - *zm); }

    T dxx() { return h2*(*xp - 2.0*(*c) + *xm); }
    T dyy() { return h2*(*yp - 2.0*(*c) + *ym); }
    T dzz() { return h2*(*zp - 2.0*(*c) + *zm); }

    T dxy() { return qrth2*(*xyp + *xym - *xm_yp - *xp_ym); }
    T dxz() { return qrth2*(*xzp + *xzm - *xm_zp - *xp_zm); }
    T dyz() { return qrth2*(*yzp + *yzm - *ym_zp - *yp_zm); }

};


#endif

