htmldir=../../html/tutorial
default: db2html
db2html:
	db2html tutorial.xml
	cp tutorial/book1.html tutorial/index.html
	cp tutorial/*html ${htmldir}
	cp `find tutorial/stylesheet-images -type f | grep -v CVS` ${htmldir}/stylesheet-images/
	cp images/*png ${htmldir}/images
	cp `find examples/born/ -type f | grep -v CVS` ${htmldir}/examples/born/
	cp `find examples/condensation/ -type f | grep -v CVS` ${htmldir}/examples/condensation/
	cp `find examples/ligand/ -type f | grep -v CVS` ${htmldir}/examples/ligand/
	cp `find examples/ion-pmf/ -type f | grep -v CVS` ${htmldir}/examples/ion-pmf/
