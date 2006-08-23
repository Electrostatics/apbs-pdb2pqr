XSLDIR=/sw/share/xml/xsl/docbook-xsl/html
XSLPATH=$(XSLDIR)/onechunk.xsl
XSLTPROCOPT=-stringparam admon.graphics 1
default: xsltproc
docbook2html:
	docbook2html --nochunks tutorial.xml
xsltproc:
	xsltproc $(XSLTPROCOPT) $(XSLPATH) tutorial.xml
