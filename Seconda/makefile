CC = g++
CFLAGS = -Wall -O3 --std=c++11



main.exe : main.o main.o
	$(CC) random.o morisetti.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
morisetti.o : morisetti.cpp morisetti.h
	$(CC) -c morisetti.cpp -o morisetti.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
