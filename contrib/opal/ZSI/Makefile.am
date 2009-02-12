## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  David Gohara
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .c .f .o

always_built_SUBDIRS = 
SUBDIRS = $(always_built_SUBDIRS) .

export PYTHONPATH:=	${SITEPKGS}

install-data-local:
	-mkdir -p ${SITEPKGS}
	${py_path} setup.py install --prefix=${prefix}
	
