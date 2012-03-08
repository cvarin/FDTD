doc_DATA = angora.pdf grid_origin.jpg random_sample.jpg column_major.jpg movie_file_layout.jpg

EXTRA_DIST = angora.pdf \
	grid_origin.eps grid_origin.jpg \
	random_sample.eps random_sample.jpg \
	column_major.eps column_major.jpg \
	movie_file_layout.eps movie_file_layout.jpg

info_TEXINFOS = angora.texi

#AM_MAKEINFOHTMLFLAGS = --no-headers --no-split
#AM_MAKEINFOHTMLFLAGS = --no-split

TEXI2HTML = texi2html
TEXI2HTML_FLAGS = --iftex --l2h --split=section
#--iftex --l2h

html: angora.html

angora.html: $(info_TEXINFOS)
	$(TEXI2HTML) $(TEXI2HTML_FLAGS) $(srcdir)/angora.texi
