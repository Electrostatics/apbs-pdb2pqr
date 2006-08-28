XSLDIR=/sw/share/xml/xsl/docbook-xsl/xhtml
XSLPATH=$(XSLDIR)/onechunk.xsl
XSLPATH=$(XSLDIR)/docbook.xsl
XSLTPROCOPT= -stringparam admon.graphics 1 \
	    -stringparam section.autolabel 1 \
	    -stringparam toc.section.depth 1 \
	    -stringparam section.label.includes.component.label 1
default: xsltproc
docbook2html:
	docbook2html --nochunks tutorial.xml
xsltproc:
	xsltproc $(XSLTPROCOPT) $(XSLPATH) tutorial.xml > index.html
