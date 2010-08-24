## ###########################################################################
## File:    apbs.spec
##
## Purpose: Spec file for building RPMS
##
## Notes:   If this is installed in the top directory the user can build a
##          full set of src and arch rpms with one command:
##
##          rpm -ta apbs-%{version}.tar.gz
##
## Author:  Stephen Bond and Nathan Baker
## ###########################################################################

Summary: Adaptive Poisson Boltzmann Solver
Name: apbs
Version: 1.3
Release: 1
License: GPL
Group: Applications/Science
Vendor: Baker Research Group, Washington University in St. Louis
Prefix: /usr/local
Buildroot: %{_topdir}/buildroot
Requires: maloc
Source: apbs-%{version}.tar.gz
URL: http://agave.wustl.edu/apbs/
Packager: Nathan Baker <baker@biochem.wustl.edu>
%description
APBS is a software package for the numerical solution of the Poisson-Boltzmann 
equation (PBE), one of the most popular continuum models for describing 
electrostatic interactions between molecular solutes in salty, aqueous media. 
Continuum electrostatics plays an important role in several areas of 
biomolecular simulation, including:

    * simulation of diffusional processes to determine ligand-protein and 
        protein-protein binding kinetics,
    * implicit solvent molecular dynamics of biomolecules ,
    * solvation and binding energy calculations to determine ligand-protein 
        and protein-protein equilibrium binding constants and aid in rational 
        drug design,
    * and biomolecular titration studies. 

APBS was designed to efficiently evaluate electrostatic properties for such 
simulations for a wide range of length scales to enable the investigation of 
molecules with tens to millions of atoms.

This software was primarily by written Nathan Baker during his graduate work 
with J. Andrew McCammon and Michael Holst and enhanced by contributions from 
several other authors. 

%package examples
Summary: Examples for APBS
Group: Applications/Science
Prefix: %{prefix}
BuildRequires: maloc
Requires: maloc, apbs = %{version}-%{release}
%description examples
This package contains examples for APBS

%package tools
Summary: Tools for APBS
Group: Applications/Science
Prefix: %{prefix}
Requires: maloc, apbs = %{version}-%{release}
%description tools
This package contains tools for APBS

%prep
%setup -n apbs-%{version}

%build
# We're assuming Intel compilers for Intel platforms.  These are the specific
# Intel platforms we'll support:
arch=%_arch
host=%_host
echo RPM VARIABLES:  ARCH: ${arch}, HOST: ${host}
%ifarch %{ix86}
  echo GOT INTEL ix86 ARCH
  %ifnarch i386 i486 i586 i686
    echo ASSUMING i786 ARCH!!!
    export CC="icc" 
    export CFLAGS="-O2 -tpp7 -static-libcxa -static -mp" 
    export CXX="icc" 
    export CXXFLAGS="-O2 -tpp7 -static-libcxa -static -mp" 
    export F77="ifort" 
    export FFLAGS="-O2 -tpp7 -static-libcxa -static -mp" 
    export LDFLAGS="-L${FETK_LIBRARY} -static-libcxa -i-static"
    ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix} --disable-shared
    make
  %else
    %ifarch i686
      echo CONFIGURING FOR i686 ARCH
      export CC="icc" 
      export CFLAGS="-O2 -tpp6 -static-libcxa -static -mp" 
      export CXX="icc" 
      export CXXFLAGS="-O2 -tpp6 -static-libcxa -static -mp" 
      export F77="ifort" 
      export FFLAGS="-O2 -tpp6 -static-libcxa -static -mp" 
      export LDFLAGS="-L${FETK_LIBRARY} -static-libcxa -i-static"
      ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix} --disable-shared
      make
    %else
      echo CONFIGURING FOR GENERIC i386 ARCH
      export CC="icc" 
      export CFLAGS="-O2 -static-libcxa -static -mp" 
      export CXX="icc" 
      export CXXFLAGS="-O2 -static-libcxa -static -mp" 
      export F77="ifort" 
      export FFLAGS="-O2 -static-libcxa -static -mp" 
      export LDFLAGS="-L${FETK_LIBRARY} -static-libcxa -i-static"
      ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix} --disable-shared
      make
    %endif
  %endif
%endif

# We're assuming Compaq compilers for Alpha/Linux platforms.  These are
# the specific Alpha platforms we'll support:
%ifarch alphaev6
  export CC='ccc'
  export CFLAGS='-O2 -arch ev6'
  export F77='fort'
  export FFLAGS='-O2 -arch ev6'
  ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix} --disable-tools
  make
%else
  %ifarch alpha
    export CC='ccc'
    export CFLAGS='-O2'
    export F77='fort'
    export FFLAGS='-O2'
    ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix} --disable-tools
    make
  %endif
%endif

# For Portland group compilers on the AMD Opteron
# NOTE: you must set the PGI_BLAS environment variable to the BLAS lib dir!

%ifarch x86_64
   export CC=pgcc
   export CFLAGS='-O2 -fastsse -fPIC'
   export F77=pgf77
   export FFLAGS='-O2 -fastsse -fPIC'
   export F77FLAGS='-O2 -fastsse -fPIC'
   export CXX=pgf77
   export CXXFLAGS='-O2 -fastsse -fPIC'
   export LDFLAGS='-Bstatic'	   
   ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix} --with-blas="-L${PGI_BLAS} -lblas"
   make
%endif

# For Itanium ia64 using intel-mkl BLAS
# NOTE: you must set the INTEL_BLAS environment variable to the BLAS lib dir!
# From intel-mkl notes: Use the following linking flag order:
#       -lmkl_lapack -lmkl_ipf -ldl

%ifarch ia64
   export CC=icc
   export CFLAGS='-fPIC  -O3 -tpp2' 
   export F77=ifort
   export FFLAGS='-fPIC -O3 -tpp2'
   export LDFLAGS='-static -static-libcxa'
   ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix} --with-blas="-L${INTEL_BLAS} -lmkl_lapack -lmkl_ipf -ldl" --with-blas-name="mkl_lapack"
   make
%endif

# For power 64
# NOTE:  There are a couple of changes you must make to compile with xlc
#        and xlf:
#        1.  In order to use static libraries you must edit the configuration
#            files for xlc (vac.cfg) and xlf (xlf.cfg).  For xlf, model the
#            stanza after xlf, but change the following libs in the
#            libraries_64 variable to the COMPLETE path to the static library:
#                 a) -lxlf90  to /path/to/libxlf90.a
#                 b) -lxlomp_ser to /path/to/libxlomp_ser.a
#                 c) -lxlfmath to /path/to/libxlfmath.a
#            and give the stanza a new name.
#        2.  Copy the xlc stanza in vac.cfg and make a new stanza, using the
#            same name you used in the xlf config file.
#        3.  Set var in FFLAGS -F:var below to this name.
#        4.  Set the BLAS environment variable and BLAS library flags to
#            point to the appropriate libraries (ATLAS or ESSL) - if using
#            ATLAS, make sure to use xlf when compiling the Fortran wrappers.

%ifarch ppc64 ppc64pseries ppc64iseries
     export CC=xlc
     export CFLAGS='-O2 -q64 -qarch=pwr4 -qtune=pwr4'
     export F77=xlf
     export FFLAGS='-q64 -qarch=pwr4 -qtune=pwr4 -Wl,-static -F:xtest'
     ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix} --with-blas="-L${BLAS} -lxlfblas -latlas" --with-blas-name="xlfblas"
     make
%endif

%install
mkdir -p ${RPM_BUILD_ROOT}/%{prefix}/apbs-%{version}
make install

cp -R examples ${RPM_BUILD_ROOT}/%{prefix}/apbs-%{version}/examples
cp -R tools  ${RPM_BUILD_ROOT}/%{prefix}/apbs-%{version}/tools

# Link the binary to the bin dir instead of the platform-specific dir
pushd ${RPM_BUILD_ROOT}/%{prefix}/bin
#ln -s */apbs apbs
popd

%clean
rm -rf ${RPM_BUILD_ROOT}

%post

%postun

%files
%defattr(-,root,root)
%{prefix}/bin
%{prefix}/lib
%{prefix}/include
%{prefix}/examples
%{prefix}/tools
%dir %{prefix}/apbs-%{version}
%doc AUTHORS COPYING INSTALL NEWS ChangeLog doc
%files examples
%defattr(-,root,root)
%{prefix}/apbs-%{version}/examples
%files tools
%defattr(-,root,root)
%{prefix}/apbs-%{version}/tools

