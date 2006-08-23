htmldir=../../html/tutorial
default: docbook2html
docbook2html:
	docbook2html --nochunks --output ${htmldir} tutorial.xml
	cp images/*png ${htmldir}/images
	cp `find examples/born/ -type f | grep -v .svn` ${htmldir}/examples/born/
	cp `find examples/condensation/ -type f | grep -v .svn` ${htmldir}/examples/condensation/
	cp `find examples/ligand/ -type f | grep -v .svn` ${htmldir}/examples/ligand/
	cp `find examples/ion-pmf/ -type f | grep -v .svn` ${htmldir}/examples/ion-pmf/
