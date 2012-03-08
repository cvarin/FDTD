EXTRA_DIST = LEGAL

dist_pkgdata_DATA = cfg/config_all.cfg \
	matlab/angora_movie.m matlab/angora_line.m matlab/angora_fieldvalue.m \
	matlab/hdf5_exists.m matlab/hdf5_read.m

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src doc

bin_PROGRAMS = angora

angora_SOURCES =
# Dummy C++ source to cause C++ linking.
nodist_EXTRA_angora_SOURCES = dummy.cxx

angora_LDFLAGS = -L$(blitzdir)/lib -L$(libconfigdir)/lib

angora_LDADD = src/libangora.la \
	-lblitz \
	-lconfig++

# needed for the HDF5 compiler
export HDF5_CXX
export HDF5_CLINKER
# The MPI wrapper is the first one around g++, so the above is not needed for MPI.
# If the MPI wrapper is no longer the innermost wrapper in the future, the environment variable
# that needs to be exported is "OMPI_CXX" (see man mpicxx)
