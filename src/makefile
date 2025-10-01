BUILD = ../build
VPATH = .:$(BUILD)
#HOME = /home/geev/
#HOME = /home/AD.NORCERESEARCH.NO/geev
HOME = /private/knss
NCFDIR = $(HOME)/prog/netcdf/netcdf-f

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .c


ifdef GFORTRAN
LD = gfortran
F90 = gfortran
F77 = gfortran
else
LD = pgf90
F90 = pgf90
F77 = pgf90
endif

ifdef GFORTRAN
# gfortran flags
FFLAGS = -c -g -mcmodel=large -O2 -std=gnu -fimplicit-none -frecursive \
          -fmax-errors=10  -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -Werror -Wunused \
          -frange-check  -pedantic -pedantic-errors -ffpe-trap=invalid,zero,overflow -Ofast
else
FFLAGS = -c -O3 -g -ffpe-trap=invalid -lineinfo \
      -Mpreprocess -Mcache_align -Mlfs -mcmodel=medium \
      -Mlarge_arrays -traceback -Ktrap=fp
endif
LINKFLAGS =

ifdef D3Q19
FFLAGS+= -DD3Q19
endif

ifdef DP
ifdef GFORTRAN
FFLAGS+= -DDOUBLE_PRECISION -fdefault-real-8 
else
FFLAGS+= -DDOUBLE_PRECISION -r8
endif
endif

ifdef CUDA
#FFLAGS+=  -gpu=lineinfo  -gpu=cc12,unroll -Minfo=accel -cuda
FFLAGS+=  -gpu=lineinfo  -gpu=cc89,unroll -Minfo=accel -cuda
#FFLAGS+=  -gpu=lineinfo  -gpu=ccall,unroll -Minfo=accel -cuda
LINKFLAGS+= -cuda
endif

ifdef MP
ifdef GFORTRAN
FFLAGS+= -fopenmp -DOPEN_MP
LINKFLAGS+= -fopenmp
else
FFLAGS+= -mp  -DOPEN_MP
LINKFLAGS+= -fopenmp
endif
endif

# netcdf
FFLAGS+= -I$(NCFDIR)/include


F77FLG =
F90FLG =

ifdef CUDA
LIBS =  -L$(HOME)/prog/nvhpc/hpc_sdk/Linux_x86_64/24.5/math_libs/lib64 -lcublas
else
LIBS =  -L$(HOME)/prog/nvhpc/hpc_sdk/Linux_x86_64/24.5/math_libs/lib64 -lnvblas
endif

ifdef DP
LIBS+= -lfftw3
else
LIBS+= -lfftw3f
endif

#netcdf
LIBS+= -L$(NCFDIR)/lib -lnetcdff

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


all: new $(TARGET) install

libs: libenkfanalysis.a -lcuda -lcublas


libenkfanalysis.a:
	cd ../../EnKF_analysis/lib; make BUILD=../../LBM/build

$(TARGET): $(INCFILES) $(OBJECTS)  timer.o
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

source:
	../bin/mksource.sh > source.files
	sed -i '/ cudafor/d' source.files
	sed -i '/ cublas/d' source.files
	sed -i '/ iso_fortran_env/d' source.files
depend:
	../bin/mkdepend.pl | sort -u > depends.file
	sed -i '/ cudafor/d' depends.file
	sed -i '/ cublas/d'  depends.file
	sed -i '/ iso_fortran_env/d' depends.file
	sed -i '/ netcdf/d' depends.file

tags: source.files
	tags.sh
	#ctags *.F90
	#f90tags.sh

include depends.file
