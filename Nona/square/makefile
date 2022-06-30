CC = g++
CFLAGS = -Wall -O3 --std=c++11

square.exe : square.o random.o genetic.o
	$(CC) random.o genetic.o square.o -o square.exe
square.o : square.cpp
	$(CC) -c square.cpp -o square.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
genetic.o : genetic.cpp genetic.h
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)

clean :
	rm *.o square.exe
