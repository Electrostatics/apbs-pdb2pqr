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
Version: 0.2.5
Release: 1
Copyright: GPL
Group: Applications/Science
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

%ifarch %{ix86}
%package tools
Summary: Tools for APBS
Group: Applications/Science
Prefix: %{prefix}
Requires: maloc, apbs = %{version}-%{release}
%description tools
This package contains tools for APBS
%endif

%prep
%setup -n apbs-%{version}

%build

# export FETK_INCLUDE=%{prefix}/lib
# export FETK_LIBRARY=%{prefix}/include
export FETK_PREFIX=`rpm -q maloc --queryformat "%{INSTALLPREFIX}"`
export FETK_INCLUDE=$FETK_PREFIX/include
export FETK_LIBRARY=$FETK_PREFIX/lib


# We're assuming Intel compilers for Intel platforms
%ifarch i686
  export CC="icc" 
  export CFLAGS="-O3 -tpp6 -static-libcxa" 
  export F77="ifc" 
  export FFLAGS="-O3 -tpp6 -static-libcxa" 
  export LDFLAGS="-static-libcxa"
  ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix}
  make
%else
  %ifarch %{ix86}
    export CC="icc" 
    export CFLAGS="-O3 -static-libcxa" 
    export F77="ifc" 
    export FFLAGS="-O3 -static-libcxa" 
    export LDFLAGS="-static-libcxa"
    ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix}
    make 
  %else
    %ifarch alpha
      export CC='ccc'
      export CFLAGS='-O2 -arch ev6'
      export F77='fort'
      export FFLAGS='-O2 -arch ev6'
      ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix} --disable-tools
      make
    %endif
  %endif
%endif

%install
mkdir -p ${RPM_BUILD_ROOT}/%{prefix}/apbs-%{version}
make install

mv examples ${RPM_BUILD_ROOT}/%{prefix}/apbs-%{version}/examples
%ifarch %{ix86}
mv tools  ${RPM_BUILD_ROOT}/%{prefix}/apbs-%{version}/tools
%endif

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
%ifarch %{ix86}
  %files tools
  %defattr(-,root,root)
  %{prefix}/apbs-%{version}/tools
%endif
