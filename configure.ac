define(LIBAPN_COPYRIGHT,[[
Copyright 2017 Jean-Pierre Flori, Jérémy Jean

This file is part of libapn.

libapn is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

libapn is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with libapn; if not, see <http://www.gnu.org/licenses/>.
]])

AC_COPYRIGHT(LIBAPN_COPYRIGHT)
AH_TOP(/*LIBAPN_COPYRIGHT*/)

AC_INIT([libapn], [1.0], [jean-pierre.flori@ssi.gouv.fr,jeremy.jean@ssi.gouv.fr])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])

dnl Automake
AM_INIT_AUTOMAKE()

dnl Libtool
LT_INIT([disable-shared])

AC_CONFIG_HEADERS(config.h:config.in)
AC_CANONICAL_HOST

${CXXFLAGS="-O3 -march=native -Wall -Wextra -pedantic"}

dnl compiler checks
AC_LANG([C++])
AC_PROG_CPP
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_14(noext)

dnl gcc builtins check
dnl AX_GCC_BUILTIN([__builtin_ffs]
dnl pthread

AC_ARG_ENABLE(pthread,
AC_HELP_STRING([--enable-pthread], [enable pthread [default=yes]]),
[case $enableval in
yes|no);;
*) AC_MSG_ERROR([bad value $enableval for --enable-pthread, need yes or no]);;
esac],
[enable_pthread=yes])

if test "$enable_pthread" = "yes"; then
AX_PTHREAD
fi

dnl debug output
AC_ARG_ENABLE(debug,
AC_HELP_STRING([--enable-debug], [enable debug [default=no]]),
[case $enableval in
yes|no);;
*) AC_MSG_ERROR([bad value $enableval for --enable-debug, need yes or no]);;
esac],
[enable_debug=no])
case $enable_debug in
yes) AC_DEFINE([DEBUG], [1], [Enable debug output]);;
no) ;;
esac

dnl Output
AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile])
AC_OUTPUT
