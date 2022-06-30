CC = g++
CFLAGS = -Wall -O3 --std=c++11

optimize.exe : optimize.o random.o
	$(CC) random.o optimize.o -o optimize.exe
optimize.o : optimize.cpp
	$(CC) -c optimize.cpp -o optimize.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o optimize.exe
