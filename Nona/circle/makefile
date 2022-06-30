CC = g++
CFLAGS = -Wall -O3 --std=c++11

circle.exe : circle.o random.o genetic.o
	$(CC) random.o genetic.o circle.o -o circle.exe
circle.o : circle.cpp
	$(CC) -c circle.cpp -o circle.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
genetic.o : genetic.cpp genetic.h
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)

clean :
	rm *.o circle.exe
