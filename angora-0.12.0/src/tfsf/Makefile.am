SUBDIRS = glquad \
	evan

AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src \
	-I$(boostdir)

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libtfsf.la

libtfsf_la_SOURCES = Ctfsf.cpp Ctfsf.h \
	Cpw.cpp Cpw.h \
	Cpw_fs.cpp Cpw_fs.h \
	Cpw_fs_inc.cpp \
	Cpw_fs_upd.cpp \
	Cpw_2l.cpp Cpw_2l.h \
	Cpw_2l_inc.cpp \
	Cpw_ml.cpp Cpw_ml.h \
	Cpw_ml_inc.cpp \
	Cpw_ml_upd.cpp \
	Cfb.cpp Cfb.h \
	Chgb.cpp Chgb.h \
	Cgsmb.cpp Cgsmb.h \
	Ckohler.cpp Ckohler.h \
	read_tfsf.cpp read_tfsf.h

# Add the libtool convenience libraries for
#  evanescent-wave and Gauss-Legendre functionalities
libtfsf_la_LIBADD = \
	evan/libevan.la \
	glquad/libglquad.la
