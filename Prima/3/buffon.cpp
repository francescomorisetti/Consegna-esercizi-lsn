#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "morisetti.h"
#include "random.h"


using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;



int n=100000;
int N=100;
int M=100;


double l=0.5;
double d=1.;

double y0; //coordinata del centro della barretta
double theta; //angolo formato dalla barretta con le ascisse uniformemente distribuito tra -1 e 1
int Nhit;
double pi[n];

double av[N];
double av2[N];

double sum_prog[N];
double su2_prog[N];
double err_prog[N];

double sumprog;

for(int j=0; j<N; j++) {
sumprog=0;

for(int k=0; k<M; k++){	//effettuo una stima di pi M volte e ne faccio la media

	Nhit=0;
	for(int i=0; i<n; i++){	//lancio la barretta 10000 volte

	y0=rnd.Rannyu();
	theta=rnd.Rannyu(0., 4*atan(1));
	
	
	if(y0<=0.5) {
	
		if(y0-deltay(l, cos(theta))<=0) Nhit++;
		
		}

	else if(y0+deltay(l, cos(theta))>=1) Nhit++;
		

				}

pi[k]=(2.*l*n)/(Nhit*d);
sumprog+=pi[k];

			}
av[j]=sumprog/N;
av2[j]=pow(av[j],2);

}


for(int i=0; i<N; i++){
	for(int j=0; j<i+1; j++){

		sum_prog[i] += av[j];
		su2_prog[i] += av2[j];
				}
	
	sum_prog[i]/=(i+1);
	su2_prog[i]/=(i+1);
	if(i==0) err_prog[i]=0;
	else err_prog[i]=sqrt((su2_prog[i] - pow(sum_prog[i],2))/i);
}

ofstream uscita;
uscita.open("pi.txt");
	for(int k=0;k<N;k++) uscita<<sum_prog[k]<<" "<<err_prog[k]<<endl;
uscita.close();



return 0;
}










