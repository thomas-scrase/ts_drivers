#Include commands common to every Makefile.am
# include $(top_srcdir)/config/makefile_templates/demo_drivers
# Name of executables that should run during the self-test
noinst_PROGRAMS= \
 one_d

# Sources for executable
one_d_SOURCES = one_d.cc


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
one_d_LDADD = -L@libdir@ -lmonodomain -lmeshes -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)


AM_CPPFLAGS += -I@includedir@