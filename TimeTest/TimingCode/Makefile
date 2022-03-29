SRCS = timer.cpp fortranTest.f90 testtimer.cpp
OBJS = timer.o fortranTest.o testtimer.o
LIBS = -lstdc++

timer: $(OBJS)
	gfortan -o timer $(OBJS) $(LIBS)

fortranTest.o: fortranTest.f90
	gfortran fortranTest.f90
testtimer.0: testtimer.cpp
	g++ testtimer.cpp
timer.o: timer.cpp
	g++ timer.cpp
clean: 
	rm -f .o .txt