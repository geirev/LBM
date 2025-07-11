BUILD = ../build
VPATH = .:$(BUILD)
HOME = /home/geev/
HOME = /home/AD.NORCERESEARCH.NO/geev

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f



LD = pgf90   #gfortran
F90 = pgf90   #gfortran
F77 = pgf90   #gfortran


#FFLAGS = -c -r8 -Mpreprocess -O2 -Mcache_align -Mlfs -mcmodel=medium

FFLAGS = -c -r8 -Mpreprocess -mp=gpu -O2 -Mcache_align -Mlfs -mcmodel=medium\
         -Minfo=accel -Mstack_arrays -Mlarge_arrays\
         -gpu=ccall,lineinfo,maxregcount:255\
         -Minfo=accel\
         -traceback -Ktrap=fp
         
ifdef CUDA
FFLAGS+= -cuda
endif

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
	@rm -f ./$(BUILD)/$*.f
	cd $(BUILD) ; $(F77) $(FFLAGS) $(F77FLG) -o $*.o ../src/$*.F



include source.files
include target.mk

INCFILES =$(INCSOURCE:.H=.h)
F90FILES =$(F90SOURCE:.F90=.f90)
F77FILES =$(F77SOURCE:.F=.f)
OBJECTS = $(F90SOURCE:.F90=.o) $(F77SOURCE:.F=.o)


all: touch libsampling.a libpreplot.a new $(TARGET) install

libs: libsampling.a libenkfanalysis.a libpreplot.a -lcuda -lcublas

libsampling.a: all.compile
	cd ../../EnKF_sampling/lib; make -f makefile.cuda  BUILD=../../LBM/build

libenkfanalysis.a:
	cd ../../EnKF_analysis/lib; make BUILD=../../LBM/build

libpreplot.a: all.compile
	cp $(HOME)/Dropbox/Progs/tecplot-binary-files/libpreplot.a $(HOME)/Dropbox/Progs/LBM/build

$(TARGET): $(INCFILES) $(OBJECTS) libsampling.a
	@mkdir -p $(BUILD)
	cd ./$(BUILD) ; $(LD) $(LINKFLAGS) -o $(TARGET) $(OBJECTS) $(LIBDIR) $(LIBS)

install:
	cp $(BUILD)/$(TARGET) $(BINDIR)

veryclean:
	cd $(BUILD) ; rm -rf *.o *.mod *.f90 *.cf90 rii_files lib*.a

new: source depend

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
