htmldir=../../html/tutorial
default: db2html
db2html:
	db2html tutorial.xml
	cp tutorial/book1.html tutorial/index.html
	cp tutorial/*html ${htmldir}
	cp tutorial/stylesheet-images/* ${htmldir}/stylesheet-images/
	cp images/*png ${htmldir}/images
	cp examples/born/* ${htmldir}/examples/born/
	cp examples/condensation/* ${htmldir}/examples/condensation/
	cp examples/ligand/* ${htmldir}/examples/ligand/
	cp examples/ion-pmf/* ${htmldir}/examples/ion-pmf/
