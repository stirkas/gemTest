#SRCS = timer.cpp fortranTest.f90 testtimer.cpp
#OBJS = timer.o fortranTest.o testtimer.o
#LIBS = -lstdc++

timer: timer.o fortranTest.o testtimer.o
	gfortan timer.o fortranTest.o testtimer.o -o timer -lstdc++

fortranTest.o: fortranTest.f90
	gfortran fortranTest.f90 -o fortranTest.o

testtimer.0: testtimer.cpp
	g++ testtimer.cpp -o testtimer.o

fortranfucntions.o: fortranfunctions.cpp
	g++ fortranfunctions.cpp -o fortranfunctions.o

timer.o: timer.cpp
	g++ timer.cpp -o timer.o

clean: 
	rm -f .o .txt