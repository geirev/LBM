BUILD = ../build
VPATH = .:$(BUILD)
HOME = /home/geev/
HOME = /home/AD.NORCERESEARCH.NO/geev

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .c



LD = pgf90   #gfortran
F90 = pgf90   #gfortran
F77 = pgf90   #gfortran


FFLAGS = -c -r8 -O2 -g -lineinfo -gpu=lineinfo \
      -Mpreprocess -Mcache_align -Mlfs -mcmodel=medium \
      -Mlarge_arrays -traceback -Ktrap=fp \
      -gpu=ccall

ifdef CUDA
FFLAGS+= -cuda
endif
ifdef MP
FFLAGS+= -mp
endif

# gfortran flags
#FFLAGS = -c -fdefault-real-8 -mcmodel=large -O2 -std=gnu -fimplicit-none -frecursive \
#          -fmax-errors=10  -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -Werror -Wunused \
#          -frange-check  -pedantic -pedantic-errors -ffpe-trap=invalid,zero,overflow -Ofast -fopenmp


F77FLG =
F90FLG =
LINKFLAGS = -fopenmp
ifdef CUDA
LINKFLAGS+= -cuda
endif

LIBS =  ./libsampling.a ./libpreplot.a  -lblas -lfftw3 -llapack

BINDIR = $(HOME)/bin


# Rules for running cpp and updating files in $(BUILD) directory
.F90.o:
	@mkdir -p $(BUILD)
	@rm -f ./$(BUILD)/$*.*
	cd $(BUILD) ; $(F90) $(FFLAGS) $(F90FLG) -o $*.o ../src/$*.F90

.F.o:
	@mkdir -p $(BUILD)
	@rm -f ./$(BUILD)/$*.*
	cd $(BUILD) ; $(F77) $(FFLAGS) $(F77FLG) -o $*.o ../src/$*.F

.c.o:
	@mkdir -p $(BUILD)
	@rm -f ./$(BUILD)/$*.*
	cd $(BUILD) ; gcc -c ../src/$*.c


include source.files
include target.mk

INCFILES =$(INCSOURCE:.H=.h)
F90FILES =$(F90SOURCE:.F90=.f90)
F77FILES =$(F77SOURCE:.F=.f)
OBJECTS = $(F90SOURCE:.F90=.o) $(F77SOURCE:.F=.o)


all: touch libsampling.a libpreplot.a new $(TARGET) install

libs: libsampling.a libenkfanalysis.a libpreplot.a -lcuda -lcublas

libsampling.a:
	cd ../../EnKF_sampling/lib; make -B -f makefile.cuda  BUILD=../../LBM/build

libenkfanalysis.a:
	cd ../../EnKF_analysis/lib; make BUILD=../../LBM/build

libpreplot.a: all.compile
	cp $(HOME)/Dropbox/Progs/tecplot-binary-files/libpreplot.a $(HOME)/Dropbox/Progs/LBM/build

$(TARGET): $(INCFILES) $(OBJECTS) libsampling.a timer.o
	@mkdir -p $(BUILD)
	cd ./$(BUILD) ; $(LD) $(LINKFLAGS) -o $(TARGET) $(OBJECTS) timer.o $(LIBDIR) $(LIBS)

install:
	cp $(BUILD)/$(TARGET) $(BINDIR)

clean:
	cd $(BUILD) ; rm -rf *.o *.mod *.f90 *.cf90 rii_files lib*.a *.bc

new: source depend

timer.o: timer.c
	cd ./$(BUILD)
	gcc -c ../src/timer.c -o ./$(BUILD)/timer.o

touch:
	touch all.compile

source:
	../bin/mksource.sh > source.files
	sed -i '/cudafor/d' source.files
	sed -i '/cublas/d' source.files

depend:
	../bin/mkdepend_linux.pl | sort -u > depends.file
	sed -i '/cudafor/d' depends.file
	sed -i '/cublas/d'  depends.file

tags: source.files
	tags.sh
	#ctags *.F90
	#f90tags.sh

include depends.file
