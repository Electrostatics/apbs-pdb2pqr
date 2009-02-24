## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  David Gohara
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .c .f .o

always_built_SUBDIRS = ZSI

SUBDIRS = $(always_built_SUBDIRS) .

install-data-local:
	-cp ${top_apbsdir}/contrib/opal/opal-py-1.9.3/wsdl/*.py ${SITEPKGS}
	-cp ${top_apbsdir}/contrib/opal/ApbsClientHelp ${SITEPKGS}
	-cp -R ${top_apbsdir}/contrib/opal/ZSI/ZSI ${SITEPKGS}

	
	
