AC_PREREQ([2.69])
AC_INIT([icon-tixi], [2.2.2], [luis.kornblueh@mpimet.mpg.de])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([src/tixi.h])

AM_INIT_AUTOMAKE([1.16.1 foreign])
dnl Enable silent rules by default with a workaround for the case when the
dnl version constraint above is relaxed:
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AM_MAINTAINER_MODE

dnl Tests for the C compiler:
dnl Fix a bug in Autoconf 2.70 (see https://savannah.gnu.org/support/?110396):
m4_if(m4_version_compare(m4_defn([AC_AUTOCONF_VERSION]), [2.70]), [0],
  [m4_pushdef([_AC_PROG_CC_C99],
        m4_bpatsubst(m4_dquote(m4_defn([_AC_PROG_CC_C99])),
          [\[ac_c_conftest_c89_program\]], [[ac_c_conftest_c99_program]]))dnl
   AC_PROG_CC
   m4_popdef([_AC_PROG_CC_C99])],
  [AC_PROG_CC])
dnl Expand the legacy macro for older version of Autoconf (2.69 and older):
m4_if(m4_version_compare(m4_defn([AC_AUTOCONF_VERSION]), [2.70]), [-1],
  [AC_PROG_CC_C99])
AS_VAR_IF([ac_cv_prog_cc_c99], [no],
  [AC_MSG_FAILURE([C compiler does not support ISO C99])])

dnl Tests for the Fortran compiler:
AC_PROG_FC
AC_LANG_PUSH([Fortran])
dnl Automake does not trace AC_FC_PP_SRCEXT/ACX_FC_PP_SRCEXT and therefore
dnl does not add their results to the compilation rules. But AC_FC_SRCEXT is
dnl traced, therefore, we temporarily redefine AC_FC_SRCEXT so it would
dnl expand to ACX_FC_PP_SRCEXT:
m4_pushdef([AC_FC_SRCEXT], [ACX_FC_PP_SRCEXT($@)])
AC_FC_SRCEXT([f90])
m4_popdef([AC_FC_SRCEXT])
dnl Handle possibly long expansions of the preprocessor macros:
ACX_FC_LINE_LENGTH([unlimited])
ACX_FC_MODULE_NAMING([AC_SUBST([FCMODEXT], ["$acx_cv_fc_module_naming_ext"])])
AM_CONDITIONAL([FCMODUC], [test x"$acx_cv_fc_module_naming_upper" = xyes])
ACX_FC_MODULE_CHECK([ISO_C_BINDING])
AC_LANG_POP([Fortran])

dnl Tests for static library generation tools:
AC_ARG_VAR([AR], [archiver command])
AC_CHECK_TOOL([AR], [ar], [:])
AS_VAR_IF([AR], [:], [AC_MSG_ERROR([the archiver command not found])])
AC_ARG_VAR([ARFLAGS], [archiver flags])
ARFLAGS=${ARFLAGS-cr}
AC_ARG_VAR([RANLIB], [archive indexer command])
AC_PROG_RANLIB

dnl Variables for pkg-config file generation:
AC_SUBST([TIXI_PKGCONF_CLIBS], ["-licon_tixi"])
AS_IF([test -n "$LDFLAGS"],
  [AS_VAR_APPEND([TIXI_PKGCONF_CLIBS], [" $LDFLAGS"])])

dnl Checks for XML library:
ACX_LANG_PACKAGE_INIT([xml2], [-I], [-L], [/include/libxml2])

tixi_save_CFLAGS=$CFLAGS
XML2_CFLAGS=${XML2_CFLAGS-$acx_c_xml2_inc_search_args}
CFLAGS="$XML2_CFLAGS $CFLAGS"
AC_CHECK_HEADER([libxml/xmlversion.h],
  [ACX_LANG_LIB_SEARCH([XML2_CLIBS], [xmlReadFile],
     [ASX_PREPEND_LDFLAGS([$acx_c_xml2_lib_search_args], [-lxml2])],
     [XML2_CLIBS=$acx_cv_c_lib_func_xmlReadFile
      AS_IF([test -n "$XML2_CLIBS"],
        [AS_VAR_APPEND([TIXI_PKGCONF_CLIBS], [" $XML2_CLIBS"])])],
     [AC_MSG_FAILURE([cannot link to libxml2 library])])],
  [AC_MSG_FAILURE([libxml/xmlversion.h not found])], [//])
CFLAGS=$tixi_save_CFLAGS

AS_IF([test -n "$LIBS"], [AS_VAR_APPEND([TIXI_PKGCONF_CLIBS], [" $LIBS"])])

AC_CONFIG_FILES([
        Makefile
        src/Makefile
        src/pkgconfig/icon-tixi.pc
])

AC_OUTPUT
