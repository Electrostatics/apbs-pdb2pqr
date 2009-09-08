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

APBS=${top_builddir}/bin/${fetk_cpu_vendor_os}/apbs

clean:
	cd 1d7h-dmso; rm -f *.out *.log io.mc
	cd 1d7i-dss; rm -f *.out *.log io.mc

test:
	./test.sh ${APBS} ${IGEN} noocd

testocd:
	./test.sh ${APBS} ${IGEN} ocd

test.serial: test


$(srcdir)/Makefile.in:

Makefile:

