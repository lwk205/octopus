## Copyright (C) 2015 D. Strubbe
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id$
##

AC_DEFUN([ACX_PATH_METIS], [

acx_external_metis=no
AC_MSG_CHECKING(for external METIS library)

# METIS is only useful in parallel.
if test x"$acx_mpi_ok" != xyes; then
  AC_MSG_RESULT([not used without MPI])
else
  acx_external_metis=no
  AC_REQUIRE([AC_PROG_CC])

  AC_ARG_WITH([metis-prefix],
    [AC_HELP_STRING([--with-metis-prefix],
    [Directory where external METIS library was installed (must be single-precision)])])

    # FIXME: accept setting by environment variables too
  case $with_metis_prefix in
    no ) acx_external_metis=disabled ;;
    "") with_metis_prefix="/usr/include" ;;
  esac

  if test x"$acx_external_metis" != xdisabled; then
  
    dnl Backup LIBS and FCFLAGS
    acx_metis_save_CFLAGS="$CFLAGS"
    acx_metis_save_LDFLAGS="$LDFLAGS"

    if test -f "$with_metis_prefix/include/metis.h"; then
      lib_path="lib"
      include_path="include"
    fi
    if test -f "$with_metis_prefix/include/metis/metis.h"; then
      lib_path="lib"
      include_path="include/metis"
    fi
    # catch bad convention in the downloadable metis version
    if test -f "$with_metis_prefix/Lib/metis.h"; then
      lib_path=""
      include_path="Lib"
    fi
    
    METIS_INCLUDE="-I$with_metis_prefix/$include_path"
    METIS_LDFLAGS="-L$with_metis_prefix/$lib_path -lmetis"

    CFLAGS="$CFLAGS $METIS_INCLUDE"
    LDFLAGS="$LDFLAGS $METIS_LDFLAGS"

    AC_LANG_SAVE
    AC_LANG_C

    AC_LINK_IFELSE([AC_LANG_PROGRAM([
#include <metis.h>
#ifdef METIS_USE_DOUBLEPRECISION
  #error METIS must be compiled in single precision for Octopus.
#endif
],[
idx_t *options;
METIS_SetDefaultOptions(options);
    ])], [acx_external_metis=yes], [])

    AC_LANG_RESTORE
    AC_MSG_RESULT([$acx_external_metis ($METIS_INCLUDE $METIS_LDFLAGS)])

    CFLAGS="$acx_metis_save_CFLAGS"
    LDFLAGS="$acx_metis_save_LDFLAGS"
  fi

  if test x"$acx_external_metis" = xno ; then
    dnl METIS was not found to link with, but is included in the distribution

    dnl We disable METIS support only if the user is requesting this explicitly
    AC_ARG_ENABLE(metis, AS_HELP_STRING([--disable-metis], [Do not compile with internal METIS domain-partitioning library.]),[acx_internal_metis=$enableval],[acx_internal_metis=yes])

    AC_MSG_CHECKING([whether METIS included in Octopus is enabled])

    AC_MSG_RESULT([$acx_internal_metis])

    if test x"$acx_internal_metis" = xyes; then
      HAVE_METIS=1
      HAVE_COMP_METIS=1
      AC_DEFINE(HAVE_METIS, 1, [This is defined when we should compile with METIS support (default).])
      AC_DEFINE(HAVE_COMP_METIS, 1, [This is defined when we link with the internal METIS library (default).])
    else
      AC_MSG_WARN(Octopus will be compiled without METIS support)
    fi

    METIS_INCLUDE=""
    LIBS_METIS_5=""
  else
    LIBS_METIS_5="$METIS_LDFLAGS"
    acx_internal_metis=no
    AC_DEFINE(HAVE_METIS,1,[This is defined when we should compile with METIS support (default).])
  fi

  AC_SUBST(METIS_INCLUDE)
  AC_SUBST(LIBS_METIS_5)
fi
])dnl ACX_PATH_METIS
