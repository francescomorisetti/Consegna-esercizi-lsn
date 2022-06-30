#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include <vector>
#include <algorithm> 



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

   

int M=100000;
int N=100;
int L=M/N;

double r[M];
int k=0;

double ave[N];
double av2[N];
double sum_prog[N];
double su2_prog[N];
double err_prog[N];

for(int i=0; i<M; i++){
	r[i]=rnd.Rannyu();
			}
double sum1=0;

for(int i=0; i<N; i++){
	sum1=0;
	for(int j=0; j<L; j++){	
		
		k=(i*L) + j;
		sum1 += r[k];
				}
	ave[i] = sum1/L;
	av2[i] = pow(ave[i],2);
			}


for(int i=0; i<N; i++){
	for(int j=0; j<i+1; j++){

		sum_prog[i] += ave[j];
		su2_prog[i] += av2[j];
				}
	
	sum_prog[i]/=(i+1);
	su2_prog[i]/=(i+1);
	if(i==0) err_prog[i]=0;
	else err_prog[i]=sqrt((su2_prog[i] - pow(sum_prog[i],2))/i);
	

			}
ofstream uscita;						
uscita.open("1.1.txt");
	for(int k=0;k<N;k++) uscita<<sum_prog[k]<<" "<<err_prog[k]<<endl;
uscita.close();	


//SECONDA PARTE

double AVE[N];
double AV2[N];
k=0.;

for(int i=0; i<M; i++){
	r[i]=rnd.Rannyu();
			}

for(int i=0; i<N; i++){
	sum1=0;
	for(int j=0; j<L; j++){	
		
		k=(i*L) + j;
		sum1 += pow((r[k]-0.5),2);
				}
	AVE[i] = sum1/L;
	AV2[i] = pow(AVE[i],2);

			}
			
double sprog[N];
double s2prog[N];
double err[N];			
			
for(int i=0; i<N; i++){
	sprog[i]=0;
	s2prog[i]=0;
	for(int j=0; j<i+1; j++){

		sprog[i] += AVE[j];
		s2prog[i] += AV2[j];
				}
	
	sprog[i]/=(i+1);
	s2prog[i]/=(i+1);
	if(i==0) err[i]=0;
	else err[i]=sqrt((s2prog[i] - pow(sprog[i],2))/i);
	

			}
						
uscita.open("1.2.txt");
	for(int k=0;k<N;k++) uscita<<sprog[k]<<" "<<err[k]<<endl;
uscita.close();

//TERZA PARTE

int a=100;
int b=10000;

double chi[a];

uscita.open("1.3.txt");

for(int i=0; i<a; i++) {
	chi[i]=0;
	double contachi[a];
	for(int k=0; k<a; k++) contachi[k]=0; 
	
	
	for(int j=0; j<b; j++) {
	
	contachi[static_cast<int>(rnd.Rannyu()*a)]++;
	
				}
				
	for(int j=0; j<a; ++j) {
	
	chi[i]+=pow( contachi[j] - b/a, 2)/(b/a); 
	
	
				}
				
	cout<<chi[i]<<endl;
	uscita<<i+1<<" "<<chi[i]<<endl;

			}

uscita.close();




   return 0;
}
