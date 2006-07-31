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

clean:
	cd actin-dimer; make clean
	cd born; make clean
	cd FKBP; make clean
	cd hca-bind; make clean
	cd ion-pmf; make clean
	cd ion-protein; make clean
	cd ionize; make clean
	cd pka-lig; make clean
	cd point-pmf; make clean
	cd protein-rna; make clean
	cd solv; make clean

test:
	cat */TESTRESULTS.log >> OLDRESULTS.log
	rm -f */TESTRESULTS.log
	cd actin-dimer; make test
	cd born; make test
	cd FKBP; make test
	cd hca-bind; make test
	cd ion-pmf; make test
	cd ion-protein; make test
	cd ionize; make test
	cd pka-lig; make test
	cd point-pmf; make test
	cd protein-rna; make test
	cd solv; make test
	cat */TESTRESULTS.log >> TESTRESULTS.log

test.serial: test

