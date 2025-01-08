include $(top_srcdir)/config/makefile_templates/doc

docfile = obtuse_corner_flow

obtuse_corner_flow.txt: intro.txt no_correction.txt sprittles.txt correction.txt outro.txt 
	@echo "Creating obtuse_corner_flow.txt"
	@cat intro.txt no_correction.txt sprittles.txt correction.txt outro.txt  > obtuse_corner_flow.txt

view: html/index.html
	open html/index.html

figures/%.gif: figures/%.pdf
	magick -density 576 $< -resize 25% $@

figures/%.eps: figures/%.pdf
	magick -density 576 $< $@
