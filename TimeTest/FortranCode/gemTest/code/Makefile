SRCS =	gem_com.f90 gem_equil.f90 gem_main.f90 gem_outd.f90 gem_fcnt.f90 gem_fft_wrapper.f90 gem_gkps_adi.f90

OBJS =	gem_com.o gem_equil.o gem_main.o gem_outd.o gem_fcnt.o gem_fft_wrapper.o gem_gkps_adi.o

LIBS = ${DFFTPACK_DIR}/libdfftpack.a /usr/lib/x86_64-linux-gnu/lapack/liblapack.so
PLIB = gem_pputil.o

F90 = mpif90 #ftn
#OPT = -FR -r8 -heap-arrays -O2 -g -traceback -check bounds
OPT = -FR -fdefault-real-8 -O3 -ffree-line-length-512 -Wall -fcheck=all -g -fbacktrace #-r8
LDFLAGS = 

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

clean:
	rm -f *.o *.lst *.mod gem_main
