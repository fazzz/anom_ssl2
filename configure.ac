#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([assl], [VERSION], [BUG-REPORT-ADDRESS])

# Checks for programs.
AC_PROG_FC([gfortran])
AC_PROG_INSTALL

# Checks for libraries.
# FIXME: Replace 'main' with a function in '-lblas':
AC_CHECK_LIB([blas], [main])
# FIXME: Replace 'main' with a function in '-llapack':
AC_CHECK_LIB([lapack], [main])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
