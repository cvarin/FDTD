
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src \
	-I$(boostdir)

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libshape.la
libshape_la_SOURCES = Cshapes.cpp Cshapes.h \
	Cshapes_excp.h \
	Cshape.cpp Cshape.h \
	Cshape_excp.h \
	read_shapes.cpp read_shapes.h
