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
Version: 0.3.0
Release: 1
Copyright: GPL
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

# export FETK_INCLUDE=%{prefix}/lib
# export FETK_LIBRARY=%{prefix}/include
export FETK_PREFIX=`rpm -q maloc --queryformat "%{INSTALLPREFIX}"`
export FETK_INCLUDE=$FETK_PREFIX/include
export FETK_LIBRARY=$FETK_PREFIX/lib

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
    export CFLAGS="-O2 -tpp7 -static-libcxa" 
    export CXX="icc" 
    export CXXFLAGS="-O2 -tpp7 -static-libcxa" 
    export F77="ifort" 
    export FFLAGS="-O2 -tpp7 -static-libcxa" 
    export LDFLAGS="-static-libcxa"
    ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix}
    make
  %else
    %ifarch i686
      echo CONFIGURING FOR i686 ARCH
      export CC="icc" 
      export CFLAGS="-O2 -tpp6 -static-libcxa" 
      export CXX="icc" 
      export CXXFLAGS="-O2 -tpp6 -static-libcxa" 
      export F77="ifort" 
      export FFLAGS="-O2 -tpp6 -static-libcxa" 
      export LDFLAGS="-static-libcxa"
      ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix}
      make
    %else
      echo CONFIGURING FOR GENERIC i386 ARCH
      export CC="icc" 
      export CFLAGS="-O2 -static-libcxa" 
      export CXX="icc" 
      export CXXFLAGS="-O2 -static-libcxa" 
      export F77="ifort" 
      export FFLAGS="-O2 -static-libcxa" 
      export LDFLAGS="-static-libcxa"
      ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix}
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


%install
mkdir -p ${RPM_BUILD_ROOT}/%{prefix}/apbs-%{version}
make install

mv examples ${RPM_BUILD_ROOT}/%{prefix}/apbs-%{version}/examples
mv tools  ${RPM_BUILD_ROOT}/%{prefix}/apbs-%{version}/tools

# Link the binary to the bin dir instead of the platform-specific dir
pushd ${RPM_BUILD_ROOT}/%{prefix}/bin
ln -s */apbs apbs
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
%dir %{prefix}/apbs-%{version}
%doc AUTHORS COPYING INSTALL NEWS ChangeLog doc
%files examples
%defattr(-,root,root)
%{prefix}/apbs-%{version}/examples
%files tools
%defattr(-,root,root)
%{prefix}/apbs-%{version}/tools
