# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

AC_INIT([minc-stuffs],[0.1],[Jason Lerch jason@phenogenomics.ca])
AC_CONFIG_SRCDIR([src/minc_displacement.c])

AM_INIT_AUTOMAKE

# The argument is the name of the generated header file.
# It is recommended to leave it as "config.h".
AC_CONFIG_HEADERS([config.h])

with_minc2=yes

# This macro eases building the package; see m4/README.
smr_WITH_BUILD_PATH


# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL

# Checks for libraries.  See m4/README.
mni_REQUIRE_BICPL

# Finish up by writing output files.
AC_CONFIG_FILES([Makefile src/Makefile perl/Makefile data/Makefile])
AC_OUTPUT
