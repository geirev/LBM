SRC_DIR := ../src
BUILD   := ../build
VPATH   := $(SRC_DIR):$(BUILD)
HOME    := $(HOME)
BINDIR  := $(HOME)/bin
NCFDIR  := $(HOME)/.netcdf/netcdf-f

DEP_FILE := depends.file
SRC_FILE := source.files

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .c


################################# COMPILER
ifdef GFORTRAN
LD  = gfortran
F90 = gfortran
F77 = gfortran
else
LD  = nvfortran
F90 = nvfortran
F77 = nvfortran
endif

ifdef MPI
LD  = mpif90
F90 = mpif90
F77 = mpif90
endif

################################# FFLAGS
ifdef GFORTRAN
FFLAGS = -c -g -mcmodel=large -O2 -std=gnu -fimplicit-none -frecursive \
          -fmax-errors=10  -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -Werror -Wunused \
          -frange-check  -pedantic -pedantic-errors -ffpe-trap=invalid,zero,overflow
FFLAGS = -c -Ofast -mcmodel=large -march=native -std=gnu -fimplicit-none -frecursive -Wno-tabs -fno-protect-parens -ffast-math
else
ifdef CUDA
FFLAGS = -c -O2 -g -lineinfo -Mpreprocess -Mcache_align -Mlfs -mcmodel=medium -Mlarge_arrays -Mfma -traceback -Ktrap=divz,inv,ovf
FFLAGS+= -gpu=lineinfo  -gpu=cc80,unroll,ptxinfo -Minfo=accel -cuda
else
FFLAGS = -c -Ofast           -Mpreprocess -Mcache_align -Mlfs -mcmodel=medium -Mlarge_arrays -Mfma -fast
endif
endif

ifdef DP
ifdef GFORTRAN
FFLAGS+= -DDOUBLE_PRECISION -fdefault-real-8
else
FFLAGS+= -DDOUBLE_PRECISION -r8
endif
endif

ifdef GFORTRAN
FFLAGS+= -DGFORTRAN
endif

ifdef MP
ifdef GFORTRAN
FFLAGS+= -fopenmp -DOPEN_MP
else
FFLAGS+= -mp  -DOPEN_MP
endif
endif

ifdef MPI
FFLAGS+= -DMPI
endif

ifdef D3Q19
FFLAGS+= -DD3Q19
endif

ifdef NETCDF
FFLAGS+= -DNETCDF
FFLAGS+= -I$(NCFDIR)/include
endif


################################# LINKFLAGS
LINKFLAGS =
ifdef CUDA
LINKFLAGS+= -cuda -mcmodel=medium -Mlarge_arrays -gpu=cc80
endif

ifdef MP
LINKFLAGS+= -fopenmp
endif

ifdef MPI
LINKFLAGS+=
endif

F77FLG =
F90FLG =

LIBS =
ifdef DP
LIBS+= -lfftw3
else
LIBS+= -lfftw3f
endif

ifdef NETCDF
LIBS+= -L$(NCFDIR)/lib -lnetcdff
endif



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


-include $(SRC_FILE)
-include target.mk

INCFILES =$(INCSOURCE:.H=.h)
F90FILES =$(F90SOURCE:.F90=.f90)
F77FILES =$(F77SOURCE:.F=.f)
OBJECTS = $(F90SOURCE:.F90=.o) $(F77SOURCE:.F=.o)

.PHONY: all install clean new source depend all touch

all: $(SRC_FILE) $(DEP_FILE) new $(TARGET) install

$(TARGET): $(INCFILES) $(OBJECTS)  timer.o
	@mkdir -p $(BUILD)
	echo
	cd ./$(BUILD) ; $(LD) $(LINKFLAGS) -o $(TARGET) $(OBJECTS) timer.o $(LIBDIR) $(LIBS)

install:
	@mkdir -p $(BINDIR)
	cp $(BUILD)/$(TARGET) $(BINDIR)

clean:
	cd $(BUILD) ; rm -rf *.o *.mod *.f90 *.cf90 rii_files lib*.a *.bc  $(SRC_FILE)

new: depend source tags

timer.o: timer.c
	cd ./$(BUILD)
	gcc -c ../src/timer.c -o ./$(BUILD)/timer.o

################################# Automatic dependency generation
# Always regenerate these before any other target
$(SRC_FILE): $(wildcard $(SRC_DIR)/*.F90) $(wildcard $(SRC_DIR)/*.F)
	@echo ">>> Updating source.files"
	../bin/mksource.sh | sort -u > $@
	sed -i '/ cudafor/d' $@
	sed -i '/ mpi.o/d' $@
ifndef NETCDF
	sed -i '/netcdf/d' $@
endif
ifndef MPI
	sed -i '/mpi/d' $@
endif

$(DEP_FILE): $(SRC_FILE)
	@echo ">>> Checking dependencies..."
	@if [ -f $@ ]; then cp $@ $@.old; fi
	../bin/mkdepend.pl | sort -u > $@
	sed -i '/ cudafor/d;/ mpi.o/d' $@
ifndef NETCDF
	sed -i '/netcdf/d' $@
endif
ifndef MPI
	sed -i '/mpi/d' $@
endif
	@if [ -f $@.old ]; then \
		if ! cmp -s $@ $@.old; then \
			echo ">>> Dependencies updated — please rerun make."; \
			rm -f $@.old; \
			false; \
		else \
			rm -f $@.old; \
		fi; \
	else \
		echo ">>> First-time dependency generation — please rerun make."; \
		false; \
	fi

stop-deps:
	$(error Stopping because depends.file was regenerated)

touch:
	touch makefile

tags: source.files
	../bin/tags.sh

-include $(DEP_FILE)
