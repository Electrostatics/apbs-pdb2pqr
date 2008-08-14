XSLDIR=/sw/share/xml/xsl/docbook-xsl/html
XSLPATH=$(XSLDIR)/onechunk.xsl
XSLTPROCOPT=-stringparam admon.graphics 1 -stringparam section.autolabel 1 -stringparam toc.section.depth 3 -stringparam section.label.includes.component.label 1
default: xsltproc
docbook2html:
	docbook2html --nochunks user-guide.xml
xsltproc:
	xsltproc $(XSLTPROCOPT) $(XSLPATH) user-guide.xml > index.html
