
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libkiss_F_fft.la

libkiss_F_fft_la_SOURCES = kiss_F_fft.c kiss_F_fft.h \
	kiss_F_fftr.c kiss_F_fftr.h\
	kiss_F_fftnd.c kiss_F_fftnd.h \
	kiss_F_fftndr.c kiss_F_fftndr.h\
	_kiss_F_fft_guts.h
