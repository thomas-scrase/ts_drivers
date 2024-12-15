#Include commands common to every Makefile.am
# include $(top_srcdir)/config/makefile_templates/demo_drivers
# Name of executables that should run during the self-test
noinst_PROGRAMS= \
 cantilever_beam \
 cube_hang \
 beam_contraction

# Sources for executable
cantilever_beam_SOURCES = cantilever_beam.cc \
                          anisotropic_test_law.h
cube_hang_SOURCES = cube_hang.cc \
                    anisotropic_test_law.h
beam_contraction_SOURCES = beam_contraction.cc \
                           anisotropic_test_law.h


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
cantilever_beam_LDADD = -L@libdir@ -lanisotropic_solid -lanisotropic_constitutive -lsolid -lconstitutive -lmeshes -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)
cube_hang_LDADD = -L@libdir@ -lanisotropic_solid -lanisotropic_constitutive -lsolid -lconstitutive -lmeshes -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)
beam_contraction_LDADD = -L@libdir@ -lanisotropic_solid -lanisotropic_constitutive -lsolid -lconstitutive -lmeshes -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)


AM_CPPFLAGS += -I@includedir@