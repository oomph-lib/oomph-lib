include $(top_srcdir)/config/makefile_templates/doc

oomph_index.txt: intro.txt oomph.idx
	cp intro.txt oomph_index.txt
	$(top_builddir)/bin/make_index.py oomph.idx  $(top_builddir) >> oomph_index.txt

docfile = oomph_index


EXTRA_DIST += intro.txt oomph.idx README_index_format
