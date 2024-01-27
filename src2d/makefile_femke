BUILD = ../build
VPATH = .:$(BUILD)

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h



LD = gfortran
CF90 = gfortran
CF77 = gfortran

FFLAGS = -c -g -fdefault-real-8 -C -u -Xlist -O2 -std=gnu -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused -Werror -O3 -fopenmp

F77FLG =
F90FLG =
LINKFLAGS = -fopenmp
CPPARCH =
CPPMODEL  =

CPPFLAGS = -traditional -P $(CPPARCH) $(CPPMODEL)

# LIBS =    ./libenkfanalysis.a ./libsampling.a -lblas -lfftw3 -llapack

LIBDIR = -L /apps/noarch/software/fftw3 -L /cm/shared/apps/blas/gcc/3.8.0 -L /cm/shared/apps/lapack/gcc/64/3.9.0
LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a /apps/noarch/software/fftw3/lib/libfftw3.a ./libenkfanalysis.a ./libsampling.a -Wl,--end-group -lpthread -lm -ldl

BINDIR = $(HOME)/bin

CPP = /usr/bin/cpp

# Rules for running cpp and updating files in $(BUILD) directory
.H.h:
	@mkdir -p $(BUILD)
	@rm -f ./$(BUILD)/$*.h
	cat MODEL.CPP $*.H | $(CPP) $(CPPFLAGS) > ./$(BUILD)/$*.h


.F90.o:
	@mkdir -p $(BUILD)
	@rm -f ./$(BUILD)/$*.f90
	@cat MODEL.CPP $*.F90 | $(CPP) $(CPPFLAGS) > ./$(BUILD)/$*.f90
	cd $(BUILD) ; $(CF90) $(FFLAGS) $(F90FLG) -o $*.o $*.f90

.F.o:
	@mkdir -p $(BUILD)
	@rm -f ./$(BUILD)/$*.f
	@cat MODEL.CPP $*.F | $(CPP) $(CPPFLAGS) > ./$(BUILD)/$*.f
	cd $(BUILD) ; $(CF77) $(FFLAGS) $(F77FLG) -o $*.o $*.f



include source.files
include target.mk

INCFILES =$(INCSOURCE:.H=.h)
F90FILES =$(F90SOURCE:.F90=.f90)
F77FILES =$(F77SOURCE:.F=.f)
OBJECTS = $(F90SOURCE:.F90=.o) $(F77SOURCE:.F=.o)


all: new $(TARGET) install

libs: libsampling.a libenkfanalysis.a

libsampling.a:
	cd ../../EnKF_sampling/lib; make BUILD=../../EnKF_MS/build

libenkfanalysis.a:
	cd ../../EnKF_analysis/lib; make BUILD=../../EnKF_MS/build

$(TARGET): $(INCFILES) $(OBJECTS)
	@mkdir -p $(BUILD)
	cd ./$(BUILD) ; $(LD) $(LINKFLAGS) -o $(TARGET) $(OBJECTS) $(LIBDIR) $(LIBS)

install:
	cp $(BUILD)/$(TARGET) $(BINDIR)

clean:
	cd $(BUILD) ; rm -rf *.o *.mod *.f90 rii_files lib*.a

new: source depend

source:
	../bin/mksource.sh > source.files

depend:
	../bin/mkdepend_linux.pl | sort -u > depends.file

tags: $(F90FILES) $(F77FILES) $(INCFILES)
	ctags *.F90
	#f90tags.sh

include depends.file
