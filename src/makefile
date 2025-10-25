BUILD = ../build
VPATH = .:$(BUILD)
#HOME = /home/geev/
HOME = /home/AD.NORCERESEARCH.NO/geev
#HOME = /private/knss
NCFDIR = $(HOME)/.netcdf/netcdf-f

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
#FFLAGS = -c -O2 -g -lineinfo -Mpreprocess -Mcache_align -Mlfs -mcmodel=medium -Mlarge_arrays -traceback
FFLAGS = -c -O2 -g -ffpe-trap=invalid -lineinfo -Mpreprocess -Mcache_align -Mlfs -mcmodel=medium -Mlarge_arrays -traceback -Ktrap=fp
endif
LINKFLAGS =

ifdef D3Q19
FFLAGS+= -DD3Q19
endif

ifdef NETCDF
FFLAGS+= -DNETCDF
endif

ifdef GFORTRAN
FFLAGS+= -DGFORTRAN
endif

ifdef DP
ifdef GFORTRAN
FFLAGS+= -DDOUBLE_PRECISION -fdefault-real-8 
else
FFLAGS+= -DDOUBLE_PRECISION -r8
endif
endif

ifdef CUDA
FFLAGS+=  -gpu=lineinfo  -gpu=cc120,unroll,ptxinfo -Minfo=accel -cuda
#FFLAGS+=  -gpu=lineinfo  -gpu=cc89,unroll -Minfo=accel -cuda
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

ifdef NETCDF
FFLAGS+= -I$(NCFDIR)/include
endif


F77FLG =
F90FLG =

ifdef CUDA
LIBS =  -L$(HOME)/prog/nvhpc/hpc_sdk/Linux_x86_64/24.5/math_libs/lib64 
else
LIBS =  -L$(HOME)/prog/nvhpc/hpc_sdk/Linux_x86_64/24.5/math_libs/lib64 
endif

ifdef DP
LIBS+= -lfftw3
else
LIBS+= -lfftw3f
endif

ifdef NETCDF
LIBS+= -L$(NCFDIR)/lib -lnetcdff
endif

BINDIR = $(HOME)/bin


# Rules for running cpp and updating files in $(BUILD) directory
.F90.o:
	@mkdir -p $(BUILD)
	@rm -f ./$(BUILD)/$*.*
	@echo
	cd $(BUILD) ; $(F90) $(FFLAGS) $(F90FLG) -o $*.o ../src/$*.F90

.F.o:
	@mkdir -p $(BUILD)
	@rm -f ./$(BUILD)/$*.*
	@echo
	cd $(BUILD) ; $(F77) $(FFLAGS) $(F77FLG) -o $*.o ../src/$*.F

.c.o:
	@mkdir -p $(BUILD)
	@rm -f ./$(BUILD)/$*.*
	@echo
	cd $(BUILD) ; gcc -c ../src/$*.c


include source.files
include target.mk

INCFILES =$(INCSOURCE:.H=.h)
F90FILES =$(F90SOURCE:.F90=.f90)
F77FILES =$(F77SOURCE:.F=.f)
OBJECTS = $(F90SOURCE:.F90=.o) $(F77SOURCE:.F=.o)


all: new $(TARGET) install

libs: libenkfanalysis.a -lcuda 


libenkfanalysis.a:
	cd ../../EnKF_analysis/lib; make BUILD=../../LBM/build

$(TARGET): $(INCFILES) $(OBJECTS)  timer.o
	@mkdir -p $(BUILD)
	echo
	cd ./$(BUILD) ; $(LD) $(LINKFLAGS) -o $(TARGET) $(OBJECTS) timer.o $(LIBDIR) $(LIBS)

install:
	cp $(BUILD)/$(TARGET) $(BINDIR)

clean:
	cd $(BUILD) ; rm -rf *.o *.mod *.f90 *.cf90 rii_files lib*.a *.bc

new: source depend tags

timer.o: timer.c
	cd ./$(BUILD)
	gcc -c ../src/timer.c -o ./$(BUILD)/timer.o

ifdef NETCDF
source:
	../bin/mksource.sh > source.files
	sed -i '/ cudafor/d' source.files
	sed -i '/ cublas/d' source.files
	sed -i '/ iso_fortran_env/d' source.files
	sed -i '/ iso_c_binding/d' source.files
else
source:
	../bin/mksource.sh > source.files
	sed -i '/ cudafor/d' source.files
	sed -i '/ cublas/d' source.files
	sed -i '/ iso_fortran_env/d' source.files
	sed -i '/m_netcdfout/d' source.files
	sed -i '/ iso_c_binding/d' source.files
endif

ifdef NETCDF
depend:
	../bin/mkdepend.pl | sort -u > depends.file
	sed -i '/ cudafor/d' depends.file
	sed -i '/ cublas/d'  depends.file
	sed -i '/ iso_fortran_env/d' depends.file
	sed -i '/ netcdf/d' depends.file
	sed -i '/ iso_c_binding/d' depends.file
else
depend:
	../bin/mkdepend.pl | sort -u > depends.file
	sed -i '/ cudafor/d' depends.file
	sed -i '/ cublas/d'  depends.file
	sed -i '/ iso_fortran_env/d' depends.file
	sed -i '/ netcdf/d' depends.file
	sed -i '/m_netcdfout/d' depends.file
	sed -i '/ iso_c_binding/d' depends.file
endif

touch:
	touch makefile

tags: source.files
	tags.sh
	#ctags *.F90
	#f90tags.sh

include depends.file
