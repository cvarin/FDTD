
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libkiss_D_fft.la

libkiss_D_fft_la_SOURCES = kiss_D_fft.c kiss_D_fft.h \
	kiss_D_fftr.c kiss_D_fftr.h\
	kiss_D_fftnd.c kiss_D_fftnd.h \
	kiss_D_fftndr.c kiss_D_fftndr.h\
	_kiss_D_fft_guts.h
