### ##########################################################################
### MALOC = < Minimal Abstraction Layer for Object-oriented C >
### Copyright (C) 1994--2000  Michael Holst
###
### This program is free software; you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by the
### Free Software Foundation; either version 2 of the License, or (at your
### option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
### See the GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License along
### with this program; if not, write to the Free Software Foundation, Inc.,
### 675 Mass Ave, Cambridge, MA 02139, USA.
###
### rcsid="$Id: maloc.spec,v 1.6 2006/02/28 22:34:53 mholst Exp $"
### ##########################################################################

## ###########################################################################
## File:    maloc.spec
##
## Purpose: Spec file for building RPMS
##
## Notes:   If this is installed in the top directory the user can build a
##          full set of src and arch rpms with one command:
##
##          rpm -ta maloc.tar.gz
##
## Author:  Stephen Bond and Michael Holst
## ###########################################################################

Summary: Minimal Abstraction Layer for Object-oriented C
Name: maloc
Version: 0.1
Release: 2
Copyright: GPL
Group: Applications/Science
Prefix: /usr/local
Buildroot: %{_topdir}/buildroot
Source: maloc-0.1-2.tar.gz
URL: http://scicomp.ucsd.edu/~mholst
Packager: Michael Holst <mholst@math.ucsd.edu>
%description
MALOC (Minimal Abstraction Layer for Object-oriented C) is a small, portable,
abstract C environment library for object-oriented C programming. MALOC is 
used as the foundation layer for a number of scientific applications, 
including MC, SG, and APBS. MALOC can be used as a small stand-alone 
abstraction environment for writing portable C programs which need access to 
resources which are typically architecture-dependent, such as INET sockets, 
timing routines, and so on. MALOC provides abstract datatypes, memory 
management routines, timing routines, machine epsilon, access to UNIX and INET 
sockets, MPI, and so on. All things that can vary from one architecture to 
another are abstracted out of an application code and placed in MALOC. To port 
the application code to a new architecture, only the small MALOC library needs 
to be ported (usually just "./configure ; make"). MALOC takes the pain of 
varying UNIX (and Win32) platforms with differing library and header 
layouts completely out of the software development picture.

%prep
%setup -n maloc

%build

%ifarch alpha
  export CC='ccc'
  export CFLAGS='-O2'
  export F77='fort'
  export FFLAGS='-O2'
  ./configure --enable-shared --prefix=${RPM_BUILD_ROOT}/%{prefix}
  make 
%endif

# For Portland group compilers on the AMD Opteron
# Note - we need to disable blas 
%ifarch x86_64 
     export CC=pgcc
     export CFLAGS='-O2 -fastsse -fPIC'
     ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix} --disable-blas
     make
%endif
     
# For Itanium ia64 
# Note - we need to disable blas
%ifarch ia64
     export CC=icc
     export CFLAGS='-O2 -fPIC'
     ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix} --disable-blas
     make
%endif

# For power 64, disabling blas
%ifarch ppc64 ppc64pseries
     export CC=xlc
     export CFLAGS="-q64 -qarch=pwr4 -qtune=pwr4"
     ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix}  --disable-blas
     make
%endif

# All others
%ifnarch alpha x86_64 ia64 ppc64 ppc64pseries
     export CC=icc
     export CFLAGS='-O2 -fPIC'
     ./configure --prefix=${RPM_BUILD_ROOT}/%{prefix}
     make 
%endif

%install
mkdir -p ${RPM_BUILD_ROOT}/%{prefix}
make install

%clean
rm -rf ${RPM_BUILD_ROOT}

%post

%postun

%files
%defattr(-,root,root)
%{prefix}/lib
%{prefix}/include
%doc AUTHORS COPYING INSTALL NEWS ChangeLog doc

