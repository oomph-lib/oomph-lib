
#To create the documentation, we need to create the html directory with its
# index file
noinst_SCRIPTS = html/index.html

# html directory depends on the doxygenified version of doc.txt
# which needs to be doxygen-ed, using the local Doxyfile
# If we have latex and dvipdf then create pdf-based documentation
#html/index.html: oomph-lib_header.html
html/index.html: raw_doxygen.h
        # Make, then kill symbolic link, but don't delete the directory itself!
#	if [ -L extra_latex_style_files ]; then rm -f extra_latex_style_files; fi
#	if [ ! -d extra_latex_style_files ]; then ln -sf $(top_srcdir)/doc/extra_latex_style_files; fi
	doxygen

EXTRA_DIST = Doxyfile


#oomph-lib_header.html:
#	$(top_builddir)/bin/build_oomph_html_header.sh $(top_builddir)

#Cleanup
clean-local:
	rm -rf *~ html latex
#	rm -rf latex html $(docfile).txt_doxygenified.h $(docfile).pdf *~ .tmp
#	rm -f oomph-lib_header.html
#	rm -f oomph-lib_footer.html
#	if [ -L extra_latex_style_files ]; then rm -f extra_latex_style_files; fi


#include $(top_srcdir)/config/makefile_templates/doc
#docfile = raw_doxygen
