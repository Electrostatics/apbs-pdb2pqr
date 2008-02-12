## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   The final makefile produced by autoconf will generate basic
##			testing results via make test.
##
##
## Author:  Todd Dolinsky
## ###########################################################################

APBS=${prefix}/bin/${fetk_cpu_vendor_os}/apbs
CLEANFILES = *.out *.log io.mc

test:
	./test.sh ${APBS}

test.serial: test


$(srcdir)/Makefile.in:

Makefile:

