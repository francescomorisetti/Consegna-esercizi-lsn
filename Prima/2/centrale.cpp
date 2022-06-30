#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include <vector>
#include <algorithm> 

/*

CC = g++
CFLAGS = -Wall -O3 --std=c++11

centrale.exe : centrale.o centrale.o
	$(CC) random.o centrale.o -o centrale.exe
centrale.o : centrale.cpp
	$(CC) -c centrale.cpp -o centrale.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o centrale.exe seed.out
*/


using namespace std;
 
int main (int argc, char *argv[]){

   Random r;
   r=r.inizializza();
   


int n=10000;
int d[4]={1,2,10,100};
double SN[n][4];
double media;

//DADO UNIFORME

for (int h=0; h<4; h++) {

	for(int i=0; i<n; i++) {
	
		media=0;
		for(int j=0; j<d[h]; j++) {
	
		media+=r.Rannyu();
		
				}
	media=media/d[h];
	SN[i][h]=media;


			}
}

ofstream uscita;						
uscita.open("unif.txt");
for (int i=0; i<n; i++) {

	for(int h=0; h<4; h++) {
	
	uscita<<SN[i][h]<<" ";

				}
uscita<<endl;
				}
				
uscita.close();				
			
//DADO ESPONENZIALE

			
for (int h=0; h<4; h++) {

	for(int i=0; i<n; i++) {
	
		media=0;
		for(int j=0; j<d[h]; j++) {
	
		media+=r.espo(1.);
		
				}
	media=media/d[h];
	SN[i][h]=media;


			}
}

						
uscita.open("exp.txt");
for (int i=0; i<n; i++) {

	for(int h=0; h<4; h++) {
	
	uscita<<SN[i][h]<<" ";

				}
uscita<<endl;
				}	
				
uscita.close();
			
				
//DADO LORENTZIANO
			
for (int h=0; h<4; h++) {

	for(int i=0; i<n; i++) {
	
		media=0;
		for(int j=0; j<d[h]; j++) {
	
		media+=r.lorentz(1.);
		
				}
	media=media/d[h];
	SN[i][h]=media;


			}
}

					
uscita.open("cauchy.txt");
for (int i=0; i<n; i++) {

	for(int h=0; h<4; h++) {
	
	uscita<<SN[i][h]<<" ";

				}
uscita<<endl;
				}				
				

uscita.close();

return 0;

}

