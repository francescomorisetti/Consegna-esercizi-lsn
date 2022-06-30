#include "morisetti.h"



using namespace std;



int main (int argc, char *argv[]){

			
int dim=100;	//dimensione random walk e numero di blocchi
int M=10000;	//numero di simulazioni
int L=M/dim;	//numero di lanci per blocco

double distanza=0;
double sum=0;
double sum2=0;
double err[dim];	
double RN[dim];
	

Random rand;
rand=rand.inizializza(rand);

	

	for (int j=0; j<dim; j++) {		

	RN[j]=0;
	sum=0;
	sum2=0;
	
	
	for (int k=0; k<dim; k++) {	
		
		distanza=0;
		
		for ( int i = 0 ; i < L ; i++) {	//calcolo la distanza per ogni step
			
		randomwalk w(dim, rand);
		distanza += moduloquadro(w.Getx(j), w.Gety(j), w.Getz(j));

						}
						
		distanza=distanza/L;
		RN[j]+=sqrt(distanza);
		
		sum2+=distanza;
		
					}
		
				
		RN[j]=RN[j]/dim;
		sum2=sum2/dim;
		
		err[j]=sqrt((sum2 - pow(RN[j],2))/(dim-1));
		
					
				}
				
ofstream uscita;

uscita.open("discreto.txt");

	uscita<<0<<" "<<0<<" "<<0<<endl;
	for(int k=0; k<dim; k++) {
	
		uscita<<k+1<<" "<<RN[k]<<" "<<err[k]<<endl;
		
}
uscita.close();	




for (int j=0; j<dim; j++) {		

	RN[j]=0;
	sum=0;
	sum2=0;
	
	
	for (int k=0; k<dim; k++) {	
		
		distanza=0;
		
		for ( int i = 0 ; i < L ; i++) {	//calcolo la distanza per ogni step
			
		continuo c(dim, rand);
		distanza += moduloquadro(c.Getx(j), c.Gety(j), c.Getz(j));

						}
						
		distanza=distanza/L;
		RN[j]+=sqrt(distanza);
		
		sum2+=distanza;
		
					}
		
				
		RN[j]=RN[j]/dim;
		sum2=sum2/dim;
		
		err[j]=sqrt((sum2 - pow(RN[j],2))/(dim-1));
		
					
				}
ofstream uscita;				
uscita.open("continuo.txt");

	uscita<<0<<" "<<0<<" "<<0<<endl;

	for(int k=0; k<dim; k++) {
	
		uscita<<k+1<<" "<<RN[k]<<" "<<err[k]<<endl;
		
}
uscita.close();	


return 0;
}
