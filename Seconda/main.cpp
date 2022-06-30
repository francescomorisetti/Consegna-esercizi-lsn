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

Coseno *mycos=new Coseno(M_PI/2.);

for(int i=0; i<M; i++){
	r[i]=rnd.Rannyu();
			}
			
double sum1=0;

for(int i=0; i<N; i++){
	sum1=0;
	for(int j=0; j<L; j++){	
		
		k=(i*L) + j;
		sum1 += (M_PI/2.)*mycos->Eval(r[k]);
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
uscita.open("1.txt");
	for(int k=0;k<N;k++) uscita<<sum_prog[k]<<" "<<err_prog[k]<<endl;

uscita.close();


//SECONDA PARTE

double AVE[N];
double AV2[N];

double sum[N];
double su2[N];
double err[N];
double R[M];


for(int i=0; i<M; i++){
	R[i]=(1.-sqrt(1.-rnd.Rannyu()));
			}
			

k=0;
for(int i=0; i<N; i++){
	sum1=0;
	for(int j=0; j<L; j++){	
		
		k=(i*L) + j;
		sum1 += (M_PI/4.)*(mycos->Eval(R[k]))/(1.-R[k]);
				}
	AVE[i] = sum1/L;
	AV2[i] = pow(AVE[i],2);
	
			}
						
for(int i=0; i<N; i++){

	sum[i]=0;
	su2[i]=0;
	for(int j=0; j<i+1; j++){

		sum[i] += AVE[j];
		su2[i] += AV2[j];
				}
	
	sum[i]/=(i+1);
	su2[i]/=(i+1);
	if(i==0) err[i]=0;
	else err[i]=sqrt((su2[i] - pow(sum[i],2))/i);
	cout<<sum[i]<<" "<<err[i]<<endl;
	
}
uscita.open("2.txt");
	for(int k=0;k<N;k++) uscita<<sum[k]<<" "<<err[k]<<endl;
uscita.close();	


return 0;
}
