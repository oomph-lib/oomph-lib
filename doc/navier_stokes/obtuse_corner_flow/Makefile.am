include $(top_srcdir)/config/makefile_templates/doc

DIRS = $(top_srcdir)/demo_drivers/navier_stokes/obtuse_corner_flow \
	   $(top_srcdir)/src/navier_stokes
	
DEPS = $(shell find $(DIRS) -name '*.h' -o -name '*.cc')

docfile = obtuse_corner_flow

figures/%.gif: figures/%.pdf
	magick -density 576 $< -resize 25% $@

figures/%.eps: figures/%.pdf
	magick -density 576 $< $@
