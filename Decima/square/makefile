MPI = mpicxx
CC = g++
CFLAGS = -Wall -O3 --std=c++11

square.exe : square.o random.o genetic.o
	$(MPI) random.o genetic.o square.o -o square.exe
square.o : square.cpp
	$(MPI) -c square.cpp -o square.o $(CFLAGS)
random.o : random.cpp random.h
	$(MPI) -c random.cpp -o random.o $(CFLAGS)
genetic.o : genetic.cpp genetic.h
	$(MPI) -c genetic.cpp -o genetic.o $(CFLAGS)

clean :
	rm *.o square.exe
