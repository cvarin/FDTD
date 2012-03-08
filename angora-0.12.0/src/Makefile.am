EXTRA_DIST = timing_chart.txt

SUBDIRS = main \
	config \
	shape \
	material \
	absorb \
	fft \
	specfun/hermite \
	specfun/gamma \
	imaging \
	nffft \
	tfsf \
	pointsources \
	recorder \
	estimator \
	waveforms

AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src/main \
	-I$(top_srcdir)/src/config \
	-I$(top_srcdir)/src/shape \
	-I$(top_srcdir)/src/material \
	-I$(top_srcdir)/src/absorb \
	-I$(top_srcdir)/src/fft \
	-I$(top_srcdir)/src/specfun/hermite \
	-I$(top_srcdir)/src/specfun/gamma \
	-I$(top_srcdir)/src/imaging \
	-I$(top_srcdir)/src/nffft \
	-I$(top_srcdir)/src/nffft/td \
	-I$(top_srcdir)/src/nffft/pd \
	-I$(top_srcdir)/src/tfsf \
	-I$(top_srcdir)/src/pointsources \
	-I$(top_srcdir)/src/recorder \
	-I$(top_srcdir)/src/estimator \
	-I$(top_srcdir)/src/waveforms \
	-I$(boostdir)

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libangora.la

libangora_la_SOURCES = angora.cpp \
	globals.cpp globals.h \
	path_utils.cpp path_utils.h \
	constants.h \
	headers.h \
	material.cpp material.h \
	material_id.h \
	angora_excp.h \
	float_comp.h \
	time_axis.cpp time_axis.h

libangora_la_LIBADD = \
	main/libmain.la \
	config/libcfgutil.la \
	shape/libshape.la \
	material/libmaterial.la \
	absorb/libabsorb.la \
	fft/libfft_d.la fft/libfft_f.la \
	specfun/hermite/libhermite.la \
	specfun/gamma/libgamma.la \
	imaging/libimg.la \
	nffft/libnffft.la \
	tfsf/libtfsf.la \
	pointsources/libpointsources.la \
	recorder/librec.la \
	estimator/libest.la \
	waveforms/libwfs.la
