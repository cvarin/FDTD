
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src \
	-I$(boostdir)

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libnffft_pd.la
libnffft_pd_la_SOURCES = Cnffft_pd.cpp  Cnffft_pd.h \
	Ctr_pd.cpp Ctr_pd.h \
	Ctr_pd_pp.cpp Ctr_pd_upd.cpp \
	Ctr_pd_fs.cpp  Ctr_pd_fs.h \
	Ctr_pd_fs_pp.cpp \
	Ctr_pd_ml.cpp  Ctr_pd_ml.h \
	Ctr_pd_ml_pp.cpp \
	Ctr_pd_ml_TL.cpp \
	Ctr_pd_2l.cpp  Ctr_pd_2l.h \
	Ctr_pd_2l_pp.cpp \
	Ctr_pd_2l_TL.cpp \
	Ctr_pd_3l.cpp  Ctr_pd_3l.h \
	Ctr_pd_3l_pp.cpp \
	Ctr_pd_3l_TL.cpp
