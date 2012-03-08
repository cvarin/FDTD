
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src \
	-I$(boostdir)

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libpointsources.la
libpointsources_la_SOURCES = Cpointsources.cpp Cpointsources.h \
	read_pointsources.cpp read_pointsources.h

