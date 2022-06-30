#include "morisetti.h"



using namespace std;



int main (int argc, char *argv[]){

double S0=100;	
double T=1.;
double K=100.;
double r=0.1;
double sigma=0.25;	

int M=10000;		//lanci totali
int N=100;		//numero blocchi
int L=M/N;	

Random rand;
rand=rand.inizializza(rand);

double S;
double call;
double put;

double sum1;
double sum2;

double avecall[N];
double av2call[N];
double aveput[N];
double av2put[N];

double sumcall_prog[N];
double su2call_prog[N];
double sumput_prog[N];
double su2put_prog[N];

double errcall_prog[N];
double errput_prog[N];

int k=0;

//PRIMA PARTE

for(int i=0; i<N; i++){
	sum1=0;
	sum2=0;
	
	for(int j=0; j<L; j++){	
		
		k=(i*L) + j;
		
		S=S0*exp((r - pow(sigma,2)/2.)*T + sigma*rand.Gauss(0.,T));
		
		call=exp(-r*T)*std::max(0., S-K);
		put=exp(-r*T)*std::max(0., K-S);
		
		sum1 +=call;
		sum2 +=put;
				}
				
	avecall[i] = sum1/L;
	av2call[i] = pow(avecall[i],2);
	
	aveput[i] = sum2/L;
	av2put[i] = pow(aveput[i],2);
			}


for(int i=0; i<N; i++){

	sumcall_prog[i]=0;
	su2call_prog[i]=0;	
	sumput_prog[i]=0;
	su2put_prog[i]=0;


	for(int j=0; j<i+1; j++){

		sumcall_prog[i] += avecall[j];
		su2call_prog[i] += av2call[j];
		
		sumput_prog[i] += aveput[j];
		su2put_prog[i] += av2put[j];
				}
	
	sumcall_prog[i]/=(i+1);
	su2call_prog[i]/=(i+1);
	
	sumput_prog[i]/=(i+1);
	su2put_prog[i]/=(i+1);
	
	
	if(i==0) { errcall_prog[i]=0; errput_prog[i]=0; }
	else { errcall_prog[i]=sqrt((su2call_prog[i] - pow(sumcall_prog[i],2))/i);
		errput_prog[i]=sqrt((su2put_prog[i] - pow(sumput_prog[i],2))/i);
		}

			}


ofstream uscita;						
uscita.open("direct.txt");
	for(int k=0;k<N;k++) {
	
		uscita<<sumcall_prog[k]<<" "<<errcall_prog[k]<<" ";
		uscita<<sumput_prog[k]<<" "<<errput_prog[k]<<endl;
		
				}
uscita.close();	


//SECONDA PARTE

double Dt=T/100.;
double t=0.;
k=0;
rand=rand.inizializza(rand);

for(int i=0; i<N; i++){
	sum1=0;
	sum2=0;
	
	for(int j=0; j<L; j++){	
		
		k=(i*L) + j;
		
		S=S0;
		
		
		for(int h=0; h<100; h++) {
		
		S=S*exp((r - pow(sigma,2)/2.)*Dt + sigma*rand.Gauss(0.,T)*sqrt(Dt));
		
		}
		
		
		call=exp(-r*T)*std::max(0., S-K);
		put=exp(-r*T)*std::max(0., K-S);
		
		sum1 +=call;
		sum2 +=put;
				}
				
	avecall[i] = sum1/L;
	av2call[i] = pow(avecall[i],2);
	
	aveput[i] = sum2/L;
	av2put[i] = pow(aveput[i],2);
			}


for(int i=0; i<N; i++){

	sumcall_prog[i]=0;
	su2call_prog[i]=0;	
	sumput_prog[i]=0;
	su2put_prog[i]=0;


	for(int j=0; j<i+1; j++){

		sumcall_prog[i] += avecall[j];
		su2call_prog[i] += av2call[j];
		
		sumput_prog[i] += aveput[j];
		su2put_prog[i] += av2put[j];
				}
	
	sumcall_prog[i]/=(i+1);
	su2call_prog[i]/=(i+1);
	
	sumput_prog[i]/=(i+1);
	su2put_prog[i]/=(i+1);
	
	
	if(i==0) { errcall_prog[i]=0; errput_prog[i]=0; }
	else { errcall_prog[i]=sqrt((su2call_prog[i] - pow(sumcall_prog[i],2))/i);
		errput_prog[i]=sqrt((su2put_prog[i] - pow(sumput_prog[i],2))/i);
		}

			}

uscita.open("sampled.txt");
	for(int k=0;k<N;k++) {
	
		uscita<<sumcall_prog[k]<<" "<<errcall_prog[k]<<" ";
		uscita<<sumput_prog[k]<<" "<<errput_prog[k]<<endl;
		
				}
uscita.close();	



return 0;
}
