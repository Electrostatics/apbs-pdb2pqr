## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   The final makefile produced by autoconf will generate tests
##			for all APBS example directories.
##
##
## Author:  Todd Dolinsky
## ###########################################################################

$(srcdir)/Makefile.in:

Makefile:

#SUBDIRS = actin-dimer alkanes born FKBP hca-bind ion-pmf ion-protein ionize pka-lig point-pmf protein-rna smpbe solv opal

clean:
	cd actin-dimer; make clean
	cd alkanes; make clean
	cd born; make clean
	cd FKBP; make clean
	cd hca-bind; make clean
	cd ion-pmf; make clean
	cd ion-protein; make clean
	cd ionize; make clean
	cd pka-lig; make clean
	cd point-pmf; make clean
	cd protein-rna; make clean
	cd smbe; make clean
	cd solv; make clean

test:
	rm -f TESTRESULTS.log
	touch TESTRESULTS.log
	cd actin-dimer; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd alkanes; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd born; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd FKBP; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd hca-bind; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd ion-pmf; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd ion-protein; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd ionize; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd pka-lig; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd point-pmf; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd protein-rna; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd smpbe; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd solv; make test; cat TESTRESULTS.log >> ../TESTRESULTS.log

testocd:
	rm -f TESTRESULTS.log
	touch TESTRESULTS.log
	cd actin-dimer; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd alkanes; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd born; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd FKBP; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd hca-bind; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd ion-pmf; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd ion-protein; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd ionize; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd pka-lig; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd point-pmf; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd protein-rna; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd smpbe; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log
	cd solv; make testocd; cat TESTRESULTS.log >> ../TESTRESULTS.log


test.serial: test

test-opal:
	cd ${prefix}/share/examples/opal; make test

