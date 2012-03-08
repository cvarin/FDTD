SUBDIRS = kiss_D_fft kiss_F_fft

AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libfft_d.la libfft_f.la

libfft_d_la_SOURCES = fft_d.cpp fft.h

libfft_f_la_SOURCES = fft_f.cpp fft.h

# Add the libtool convenience libraries for
#  KissFFT
libfft_d_la_LIBADD = kiss_D_fft/libkiss_D_fft.la
libfft_f_la_LIBADD = kiss_F_fft/libkiss_F_fft.la
