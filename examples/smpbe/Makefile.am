## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   The final makefile produced by autoconf will generate basic
##			testing results via make test.
##
##
## Author:  David Gohara
## ###########################################################################

APBS=${top_builddir}/bin/${fetk_cpu_vendor_os}/apbs
IGEN=${top_builddir}/tools/manip/inputgen.py
CLEANFILES = *-PE* *.out *.log io.mc*

test:
	./test.sh ${APBS} ${IGEN}

test.serial: test


$(srcdir)/Makefile.in:

Makefile:

