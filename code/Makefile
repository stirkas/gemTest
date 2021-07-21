SRCS =	gem_com.f90 gem_equil.f90 gem_main.f90 gem_outd.f90 gem_fcnt.f90 gem_fft_wrapper.f90 gem_gkps_adi.f90

OBJS =	gem_com.o gem_equil.o gem_com_c.o gem_equil_c.o spec_c.o gem_main.o gem_outd.o gem_fcnt_c.o gem_fcnt_f.o gem_fft_wrapper.o gem_gkps_adi.o
LSTDC++ = /usr/lib/gcc/x86_64-linux-gnu/7/libstdc++.so#this is trying to link the c++ library 
LIBS = /usr/local/lib/libdfftpack.a /usr/lib/x86_64-linux-gnu/liblapack.so -lstdc++
PLIB = gem_pputil.o

F90 = mpif90 #ftn
C++ = g++
#OPT = -FR -r8 -heap-arrays -O2 -g -traceback -check bounds
OPT = -FR -fdefault-real-8 -O2 -ffree-line-length-512 -fcheck=all -g -fbacktrace #-Wall
OPT_C = -O2 -std=gnu++17 #need to make sure that this works for c++14 on the cray computers
LDFLAGS = 

#all : gem

gem_main: ${OBJS}
	$(F90)  -o gem_main $(OPT) $(OBJS) $(PLIB) $(LIBS) 

gem_pputil.o: gem_pputil.f90
	$(F90) -c $(OPT) gem_pputil.f90

gem_com.o: gem_com.f90 gem_pputil.o
	$(F90) -c $(OPT) gem_com.f90

gem_com_c.o: gem_com_c.cpp
	$(C++) -c $(OPT_C) gem_com_c.cpp

gem_equil.o: gem_equil.f90 gem_pputil.o
	$(F90) -c $(OPT) gem_equil.f90

gem_equil_c.o: gem_equil_c.cpp
	$(C++) -c $(OPT_C) gem_equil_c.cpp

gem_gkps_adi.o: gem_gkps_adi.f90 gem_com.f90 gem_equil.f90 gem_pputil.f90
	$(F90) -c $(OPT) gem_gkps_adi.f90

gem_main.o: gem_main.f90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o gem_gkps_adi.o gem_fcnt_c.o
	$(F90) -c $(OPT) gem_main.f90

gem_outd.o: gem_outd.f90 gem_fft_wrapper.o gem_pputil.o gem_com.o gem_equil.o 
	$(F90) -c $(OPT) gem_outd.f90

#gem_fcnt.o: gem_fcnt.f90
#	$(F90) -c $(OPT) gem_fcnt.f90

gem_fcnt_c.o: gem_fcnt_c.cpp#this is trying to compile the cpp file
	$(C++) -c $(OPT_C) gem_fcnt_c.cpp

gem_fcnt_f.o: gem_fcnt_f.f90#this is tryint to compile the f90 file
	$(F90) -c $(OPT) gem_fcnt_f.f90

gem_fft_wrapper.o: gem_fft_wrapper.f90
	$(F90) -c $(OPT) gem_fft_wrapper.f90

spec_c.o: spec_c.cpp
	$(C++) -c $(OPT_C) spec_c.cpp

clean:
	rm -f *.o *.lst *.mod gem_main
