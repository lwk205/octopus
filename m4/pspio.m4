## Copyright (C) 2012 M. Oliveira
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

AC_DEFUN([ACX_PSPIO], [
acx_pspio_ok=no

dnl Backup LIBS and FCFLAGS
acx_pspio_save_LIBS="$LIBS"
acx_pspio_save_FCFLAGS="$FCFLAGS"

if test "$PKGCONFIG" != ""; then
  PSPIO_PREFIX=`$PKGCONFIG --variable=prefix pspio`
else
  PSPIO_PREFIX=/usr
fi
dnl Check if the library was given in the command line
AC_ARG_WITH(pspio-prefix, [AS_HELP_STRING([--with-pspio-prefix=DIR], [Directory where pspio was installed.])],[],[with_pspio_prefix=$PSPIO_PREFIX])
case $with_pspio_prefix in
  no ) acx_pspio_ok=disable ;;
  *) LIBS_PSPIO="-L$with_pspio_prefix/lib -lfpspio -lpspio"; FCFLAGS_PSPIO="$ax_cv_f90_modflag$with_pspio_prefix/include" ;;
esac

testprog="AC_LANG_PROGRAM([],[
    use fpspio_m

    type(fpspio_pspdata_t) :: pspdata
    call fpspio_pspdata_free(pspdata)])"

dnl The tests
if test "$acx_pspio_ok" = no; then
  AC_MSG_CHECKING([for pspio])
  # If the location has been passed with --with-pspio-prefix just test this
  if test "$LIBS_PSPIO"; then
    pspio_fcflags="$FCFLAGS_PSPIO"; pspio_libs="$LIBS_PSPIO"
    FCFLAGS="$pspio_fcflags $acx_pspio_save_FCFLAGS $GSL_CFLAGS"
    LIBS="$pspio_libs $acx_pspio_save_LIBS $GSL_LIBS"
    AC_LINK_IFELSE($testprog, [acx_pspio_ok=yes; FCFLAGS_PSPIO="$pspio_fcflags"; LIBS_PSPIO="$pspio_libs"], [])
  else
    pspio_libs="-lpspio_fortran -lpspio"
    FCFLAGS="$pspio_fcflags $acx_pspio_save_FCFLAGS $GSL_CFLAGS"
    LIBS=" $acx_pspio_save_LIBS $pspio_libs $GSL_LIBS"
    AC_LINK_IFELSE($testprog, [acx_pspio_ok=yes; FCFLAGS_PSPIO="$pspio_fcflags"; LIBS_PSPIO="$pspio_libs"], [])
  fi
  AC_MSG_RESULT([$acx_pspio_ok ($FCFLAGS_PSPIO $LIBS_PSPIO)])
fi

AC_SUBST(FCFLAGS_PSPIO)
AC_SUBST(LIBS_PSPIO)
FCFLAGS="$acx_pspio_save_FCFLAGS"
LIBS="$acx_pspio_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_pspio_ok" = xyes; then
  AC_DEFINE(HAVE_PSPIO,1,[Defined if you have the PSPIO library.])
else
  LIBS_PSPIO=""
  FCFLAGS_PSPIO=""
fi
])dnl ACX_PSPIO
