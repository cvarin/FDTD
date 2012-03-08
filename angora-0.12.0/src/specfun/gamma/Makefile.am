
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libgamma.la
libgamma_la_SOURCES = gamma.c gamma.h \
	polevl.c const.c mtherr.c isnan.c \
	mconf.h
