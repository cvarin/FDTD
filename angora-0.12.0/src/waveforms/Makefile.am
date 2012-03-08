
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src \
	-I$(boostdir)

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libwfs.la
libwfs_la_SOURCES = Cwfs.cpp Cwfs.h \
	Cwfs_excp.h \
	Cwf.h \
	Cwf_shared_ptr.h \
	Cwf_gauss.cpp Cwf_gauss.h \
	read_wfs.cpp read_wfs.h
