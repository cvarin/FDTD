AC_INIT([angora],[0.12.0],[bugs@angorafdtd.org])

AC_CONFIG_SRCDIR([src/constants.h])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_AUX_DIR([.])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE

# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
LT_INIT

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h sys/timeb.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor mkdir pow sqrt])

# check if MPI parallel I/O should be used
AC_ARG_ENABLE([mpi-io],
  [AS_HELP_STRING([--enable-mpi-io],
   [enables the MPI parallel I/O support [default=disabled]])],
  [AC_DEFINE([MPI_USE_PARALLEL_IO], [], [Enables MPI I/O support.])],
  [])

# location of the MPI installation directory
dnl Define a macro for introducing the MPI preprocessor variable
m4_define([ANGORA_DEFN_MPI_DISABLE],[AC_DEFINE([MPI_DISABLE], [], [Disables the message-passing-interface (MPI) support.])])
AC_ARG_WITH([mpi],
  [AS_HELP_STRING([--with-mpi@<:@=DIR@:>@],
   [enable MPI support (with binary directory DIR, if not in path)])],
  [if test "$with_mpi" = "yes" ; then
       AC_SUBST([CXX],[mpicxx])
     else
       if test "$with_mpi" != "no" ; then
           AC_SUBST([CXX],[$with_mpi/bin/mpicxx])
	   else
           ANGORA_DEFN_MPI_DISABLE
       fi
     fi],
  [ANGORA_DEFN_MPI_DISABLE])

# location of the HDF5 installation directory
#C++ compiler and linker under the h5c++ wrapper are what's been already set (g++, mpicxx, etc)
AC_SUBST([HDF5_CXX],[$CXX])dnl
AC_SUBST([HDF5_CLINKER],[$CXX])dnl
# Save the current value of CXX, to be restored later if HDF5 is disabled
CXX_SAVEFIRST=$CXX
AC_SUBST([CXX],[h5c++])dnl
AC_ARG_WITH([hdf5],
  [AS_HELP_STRING([--with-hdf5@<:@=DIR@:>@],
   [enable HDF5 support (with binary directory DIR, if not in path)])],
  [
  if test "$with_hdf5" != "no" ; then
    if test "$with_hdf5" != "yes" ; then
      AC_SUBST([CXX],[$with_hdf5/bin/h5c++])dnl
	fi
  else
    AC_DEFINE([HDF5_DISABLE], [], [Disables the HDF5 file format support.])dnl Define preprocessor variable for disabling HDF5
    AC_SUBST([CXX],[$CXX_SAVEFIRST])dnl Restore the CXX variable
  fi],
  [])dnl HDF5 enabled by default

# location of the blitz++ installation directory
AC_ARG_WITH([blitz],
  [AS_HELP_STRING([--with-blitz@<:@=DIR@:>@],
   [location of blitz++ root directory, if not in path])],
  [if test "$with_blitz" != "yes" ; then
         if test "$with_blitz" != "no" ; then
             AC_SUBST([blitzdir],[$with_blitz])
         else
             AC_MSG_ERROR([blitz++ functionality cannot be disabled])
         fi
     fi],
  [])

# location of the libconfig installation directory
AC_ARG_WITH([config],
  [AS_HELP_STRING([--with-config@<:@=DIR@:>@],
   [location of libconfig++ root directory, if not in path])],
  [if test "$with_config" != "yes" ; then
         if test "$with_config" != "no" ; then
             AC_SUBST([libconfigdir],[$with_config])
         else
             AC_MSG_ERROR([libconfig++ functionality cannot be disabled])
         fi
     fi],
  [])

# location of Boost directory
AC_ARG_WITH([boost],
  [AS_HELP_STRING([--with-boost@<:@=DIR@:>@],
   [location of Boost root directory, if not in path])],
  [if test "$with_boost" != "yes" ; then
         if test "$with_boost" != "no" ; then
             AC_SUBST([boostdir],[$with_boost])
         else
             AC_MSG_ERROR([Boost functionality cannot be disabled])
         fi
     fi],
  [])

# location of the Angora template config file
dnl See the autoconf documentation at http://www.gnu.org/software/autoconf/manual/autoconf.html#Preset-Output-Variables, section 4.8.1
dnl "...coupled with the fact that configure is always run from the top build directory,
dnl  it is sufficient to use just ‘$srcdir’ instead of ‘$top_srcdir’."
AC_DEFINE_UNQUOTED([angora_cfgdir], ["$srcdir/cfg"], [Location of the Angora template config file.])

AC_CONFIG_FILES(Makefile \
  doc/Makefile \
  src/Makefile \
  src/main/Makefile \
  src/config/Makefile \
  src/shape/Makefile \
  src/material/Makefile \
  src/absorb/Makefile \
  src/fft/Makefile \
  src/fft/kiss_D_fft/Makefile \
  src/fft/kiss_F_fft/Makefile \
  src/specfun/hermite/Makefile \
  src/specfun/gamma/Makefile \
  src/pointsources/Makefile \
  src/waveforms/Makefile \
  src/recorder/Makefile \
  src/estimator/Makefile \
  src/imaging/Makefile \
  src/nffft/Makefile \
  src/nffft/td/Makefile \
  src/nffft/pd/Makefile \
  src/tfsf/Makefile \
  src/tfsf/glquad/Makefile \
  src/tfsf/evan/Makefile)

AC_OUTPUT
