SRCS =	gem_com.f90 gem_equil.f90 gem_main.f90 gem_outd.f90 gem_fcnt.f90 gem_fft_wrapper.f90 gem_gkps_adi.f90

OBJS =	gem_com.o gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_fft_wrapper.o gem_gkps_adi.o

PLIB = gem_pputil.o
INCS = -I"${OMP_INC_DIR}"

#Set debug flag to 0 so it has to be forced to 1 to get debug builds.
DEBUG = 0

#Get the location of the dfftpack library using git root dir.
ROOT_DIR     = $(shell git rev-parse --show-toplevel)
DFFTPACK_DIR = $(ROOT_DIR)/dfftpack

#Get a list of loaded modules from Cori.
ifeq (${NERSC_HOST}, cori)
   MODULE_LIST = ${LOADEDMODULES}
endif

#Figure out which compiler to used based on loaded module list. Default to GNU.
CRAY_MODULE      = PrgEnv-cray
INTEL_MODULE     = PrgEnv-intel
INTEL_MODULE_MPI = impi
#CRAY
ifeq ($(findstring $(CRAY_MODULE), $(MODULE_LIST)), $(CRAY_MODULE))
    F90 = ftn
    C++ = CC

    OPT_COMMON   = -f free -s real64 -h omp
    OPT_C_COMMON = -fopenmp 
    ifeq ($(DEBUG),0)
        OPT   = $(OPT_COMMON)   -O3
        OPT_C = $(OPT_C_COMMON) -O3 
    else
        OPT_DEBUG_COMMON = -O2 -g
        OPT   = $(OPT_COMMON)   $(OPT_DEBUG_COMMON) -O negmsgs
        OPT_C = $(OPT_C_COMMON) $(OPT_DEBUG_COMMON)
    endif

    LIBS += $(DFFTPACK_DIR)/cray/libdfftpack.a
    LDFLAGS = 

#INTEL
else ifeq ($(findstring $(INTEL_MODULE), $(MODULE_LIST)), $(INTEL_MODULE))
    ifeq ($(findstring $(INTEL_MPI_MODULE), $(MODULE_LIST)), $(INTEL_MPI_MODULE)) #Nested if since make doesn't support AND conditionals...
        F90 = mpiifort
        C++ = mpiicpc

        OPT_COMMON   = -FR -r8 -mkl -qopenmp
        OPT_C_COMMON = -mkl -qopenmp
        ifeq ($(DEBUG),0)
            OPT   = $(OPT_COMMON)   -O3
            OPT_C = $(OPT_C_COMMON) -O3
        else
            OPT_DEBUG_COMMON = -O2 -g -heap-arrays -traceback
            OPT   = $(OPT_COMMON)   $(OPT_DEBUG_COMMON) -check bounds
            OPT_C = $(OPT_C_COMMON) $(OPT_DEBUG_COMMON)
        endif

        LIBS += $(DFFTPACK_DIR)/intel/libdfftpack.a
        LDFLAGS = 
    endif

#GNU
else
    F90 = mpif90
    C++ = mpic++

    OPT_COMMON   = -FR -fdefault-real-8 -ffree-line-length-512
    OPT_C_COMMON = $(INCS)
    ifeq ($(DEBUG),0)
        OPT   = $(OPT_COMMON)   -O3
        OPT_C = $(OPT_C_COMMON) -O3
    else
        OPT_DEBUG_COMMON = -O0 -g -fcheck=all -fbacktrace -Wall
        OPT   = $(OPT_COMMON)   $(OPT_DEBUG_COMMON)
        OPT_C = $(OPT_C_COMMON) $(OPT_DEBUG_COMMON)
    endif

    LIBS += $(DFFTPACK_DIR)/gnu/libdfftpack.a -llapack
    LDFLAGS = 
endif

#all : gem

gem_main: gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_pputil.o gem_com.o gem_fft_wrapper.o gem_gkps_adi.o
	$(F90)  -o gem_main $(OPT) $(OBJS) $(PLIB) $(LIBS) 

gem_pputil.o: gem_pputil.f90
	$(F90) -c $(OPT) gem_pputil.f90

gem_com.o: gem_com.f90 gem_pputil.o
	$(F90) -c $(OPT) gem_com.f90

gem_equil.o: gem_equil.f90 gem_pputil.o
	$(F90) -c $(OPT) gem_equil.f90

gem_gkps_adi.o: gem_gkps_adi.f90 gem_com.f90 gem_equil.f90 gem_pputil.f90
	$(F90) -c $(OPT) gem_gkps_adi.f90

gem_main.o: gem_main.f90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o gem_gkps_adi.o
	$(F90) -c $(OPT) gem_main.f90

gem_outd.o: gem_outd.f90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o
	$(F90) -c $(OPT) gem_outd.f90

gem_fcnt.o: gem_fcnt.f90
	$(F90) -c $(OPT) gem_fcnt.f90

gem_fft_wrapper.o: gem_fft_wrapper.f90
	$(F90) -c $(OPT) gem_fft_wrapper.f90
	
debug:
	make DEBUG=1

install:
	cp $(ROOT_DIR)/code/gem_main $(ROOT_DIR)/testcase

clean:
	rm -f *.o *.lst *.mod gem_main
