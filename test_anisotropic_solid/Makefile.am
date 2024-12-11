#Include commands common to every Makefile.am
# include $(top_srcdir)/config/makefile_templates/demo_drivers
# Name of executables that should run during the self-test
noinst_PROGRAMS= \
 cubic

# Sources for executable
cubic_SOURCES = cubic.cc \
                anisotropic_test_law.h


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
cubic_LDADD = -L@libdir@ -lanisotropic_solid -lanisotropic_constitutive -lsolid -lconstitutive -lmeshes -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)


AM_CPPFLAGS += -I@includedir@