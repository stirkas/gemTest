#SRCS = timer.cpp fortranTest.f90 testtimer.cpp
#OBJS = timer.o fortranTest.o testtimer.o
#LIBS = -lstdc++

timer: timer.o fortranfunction.o testtimer.o fortranTest.o 
	g++ -o timer timer.o fortranfunction.o testtimer.o fortranTest.o -lstdc++ 

fortranTest.o: fortranTest.f90
	gfortran -c fortranTest.f90

testtimer.o: testtimer.cpp
	g++ -c testtimer.cpp

fortranfucntion.o: fortranfunction.cpp
	g++ -c fortranfunctions.cpp

timer.o: timer.cpp
	g++ -c timer.cpp  

clean: 
	rm -f .o .txt