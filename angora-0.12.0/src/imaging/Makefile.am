
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src \
	-I$(boostdir)

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libimg.la

libimg_la_SOURCES = Cimgs.cpp Cimgs.h \
	Cimg.cpp Cimg.h \
	Cimg_pp.cpp \
	read_imaging.cpp read_imaging.h
