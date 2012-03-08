SUBDIRS = td \
	pd

AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src \
	-I$(boostdir)

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libnffft.la

libnffft_la_SOURCES = read_nffft_td.cpp read_nffft_td.h \
	read_nffft_pd.cpp read_nffft_pd.h

# Add the libtool convenience libraries for
#  time-domain and phasor-domain NFFFT
libnffft_la_LIBADD = \
	td/libnffft_td.la \
	pd/libnffft_pd.la
