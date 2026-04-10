SRC_DIR := ../src
BUILD   := ../build
VPATH   := $(SRC_DIR):$(BUILD)
HOME    := $(HOME)
BINDIR  := $(HOME)/bin
NCFDIR  := $(HOME)

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



##################################################################
# 1. AUTOMATIC DETECTION of system and architecture
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# If we are on macOS, we MUST use GFortran options, not NVIDIA
ifeq ($(UNAME_S), Darwin)
    GFORTRAN = 1
endif

##################################################################
# 2. COMPILER FLAG LOGIC
ifdef GFORTRAN
    # Base flags supported by all gfortran versions
    # Note: We use -cpp instead of -Mpreprocess for gfortran
    BASE_FLAGS = -c -Ofast -std=gnu -fimplicit-none -frecursive -Wno-tabs -fno-protect-parens -ffast-math -cpp

    ifeq ($(UNAME_S), Darwin)
        # --- MACOS SPECIFIC ---

        # 1. Fix Architecture Flag
        # Apple Silicon (arm64) needs -mcpu=native. Intel needs -march=native.
        ifeq ($(UNAME_M), arm64)
             ARCH_FLAG = -mcpu=native
        else
             ARCH_FLAG = -march=native
        endif

        # 3. Construct Final Flags
        # -fallow-argument-mismatch is CRITICAL for gfortran >= 10
        FFLAGS = $(BASE_FLAGS) $(ARCH_FLAG) -fallow-argument-mismatch -mtune=native #$(NC_FLAGS)
    else
        # --- LINUX SPECIFIC ---
        FFLAGS = $(BASE_FLAGS) -march=native -mcmodel=large
    endif

else
    # --- nvfortran and CUDA flags
    ifdef CUDA
        FFLAGS = -c -O2 -g -lineinfo -Mpreprocess -Mcache_align -Mlfs -mcmodel=medium \
                 -Mlarge_arrays -Mfma -traceback -Ktrap=divz,inv,ovf \
                 -gpu=lineinfo -gpu=cc120,unroll,ptxinfo -Minfo=accel -cuda
    else
        FFLAGS = -c -Ofast -Mpreprocess -Mcache_align -Mlfs -mcmodel=medium -Mlarge_arrays -Mfma -fast
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
LINKFLAGS+= -cuda -mcmodel=medium -Mlarge_arrays -gpu=cc120
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
	../bin/mksource.sh  > $@
	sed -i.bak '/ cudafor/d' $@  && rm $@.bak
	sed -i.bak '/ mpi.o/d' $@ && rm $@.bak
	sed -i.bak '/netcdf.o/d' $@ && rm $@.bak
ifndef MPI
	sed -i.bak '/mpi/d' $@ && rm $@.bak
endif

$(DEP_FILE): $(SRC_FILE)
	@echo ">>> Checking dependencies..."
	@if [ -f $@ ]; then cp $@ $@.old; fi
	../bin/mkdepend.pl | sort -u > $@
	sed -i.bak '/ cudafor/d;/ mpi.o/d' $@ && rm $@.bak
	sed -i.bak '/netcdf.o/d' $@ && rm $@.bak
ifndef MPI
	sed -i.bak '/mpi/d' $@ && rm $@.bak
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
