
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src \
	-I$(boostdir)

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libnffft_td.la
libnffft_td_la_SOURCES = Cnffft_td.cpp  Cnffft_td.h \
	Ctr_td.cpp Ctr_td.h \
	Ctr_td_fs.cpp Ctr_td_fs.h \
	Ctr_td_fs_pp.cpp \
	Ctr_td_fs_upd.cpp \
	Ctr_td_3l.cpp Ctr_td_3l.h \
	Ctr_td_3l_pp.cpp \
	Ctr_td_3l_upd.cpp \
	Ctr_td_3l_TL.cpp
