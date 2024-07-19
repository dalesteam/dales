lib_LIBRARIES = libicon_tixi.a

include_HEADERS =        \
        tixi.h           \
        tixiData.h       \
        tixiInternal.h   \
        tixiUtils.h      \
        tixi_version.h   \
        uidHelper.h      \
        xpathFunctions.h

if FCMODUC
tixi_mod = TIXI.$(FCMODEXT)
else !FCMODUC
tixi_mod = tixi.$(FCMODEXT)
endif !FCMODUC

nodist_include_HEADERS = $(tixi_mod)

pkgconfigdir = $(libdir)/pkgconfig
dist_noinst_DATA = pkgconfig/icon-tixi.pc.in
nodist_pkgconfig_DATA = pkgconfig/icon-tixi.pc

libicon_tixi_a_SOURCES = \
        tixi.f90         \
        tixi.h           \
        tixiData.h       \
        tixiImpl.c       \
        tixiInternal.c   \
        tixiInternal.h   \
        tixiUtils.c      \
        tixiUtils.h      \
        tixi_version.h   \
        uidHelper.c      \
        uidHelper.h      \
        xpathFunctions.c \
        xpathFunctions.h

AM_CFLAGS = $(XML2_CFLAGS)

$(tixi_mod): tixi.$(OBJEXT)
	@if test ! -f $@; then rm -f $<; $(MAKE) $<; fi
MOSTLYCLEANFILES = $(tixi_mod)

clean-local:
	-rm -f *.i

