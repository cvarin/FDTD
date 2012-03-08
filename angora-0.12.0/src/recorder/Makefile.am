
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src \
	-I$(boostdir)

AM_LDFLAGS = -static

noinst_LTLIBRARIES = librec.la

librec_la_SOURCES = Crecorder.cpp Crecorder.h \
	Crec.h \
	Crecorder_movie.cpp Crecorder_movie.h \
	Crecorder_line.cpp Crecorder_line.h \
	Crecorder_fieldvalue.cpp Crecorder_fieldvalue.h \
	read_recorder.cpp read_recorder.h

