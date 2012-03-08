
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libmain.la
libmain_la_SOURCES = init.cpp init.h \
	initgeom.cpp initgeom.h \
	place_obj.cpp place_obj.h \
	placegeom.cpp placegeom.h \
	matfile.h \
	random.h \
	parallel.cpp parallel.h \
	update.cpp update.h \
	read_geom.cpp read_geom.h \
	read_basedirs.cpp read_basedirs.h \
	read_multgridinfo.cpp read_multgridinfo.h \
	read_global.cpp read_global.h
