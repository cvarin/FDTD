
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libmaterial.la
libmaterial_la_SOURCES = Cmats.cpp Cmats.h \
	Cmat_excp.h \
	read_mat.cpp read_mat.h
