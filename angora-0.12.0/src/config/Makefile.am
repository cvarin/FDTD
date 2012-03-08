
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(libconfigdir)/include \
	-I$(top_srcdir)/src

AM_LDFLAGS = -static

#binary that reads the reference config file and creates the source file for the built-in config object
noinst_PROGRAMS = createcfg
createcfg_SOURCES = createcfg.cpp
createcfg_LDFLAGS = -L$(libconfigdir)/lib
createcfg_LDADD = -lconfig++

# libvalset.la is the convenience library holding the valid Angora settings
noinst_LTLIBRARIES = libcfgutil.la
libcfgutil_la_SOURCES = cfg_utils.cpp cfg_utils.h valid_settings.cpp

valid_settings.cpp: valid_settings.cc

# built source that creates the config object holding the valid Angora settings
BUILT_SOURCES = valid_settings.cc
CLEANFILES = valid_settings.cc
valid_settings.cc: createcfg createcfg.cpp $(top_srcdir)/cfg/config_all.cfg
	./createcfg $(top_srcdir)/cfg/config_all.cfg
