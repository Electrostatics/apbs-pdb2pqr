htmldir=../../html/user-guide
default: db2html
db2html:
	db2html user-guide.xml
	cp user-guide/book1.html user-guide/index.html
	mv user-guide/*.html ${htmldir}
	mv user-guide/stylesheet-images/* ${htmldir}/stylesheet-images/
	cp -rf images/*.gif ${htmldir}/images
xmlto:
	xmlto -o ${htmldir} html user-guide.xml
	cp -rf images/*.gif ${htmldir}/images

