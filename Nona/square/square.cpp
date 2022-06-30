#include "genetic.h"

int main() {

Random r;
r=r.inizializza(r);

genetic g(r);		//Crea la popolazione iniziale e le posizioni spaziali
g.Fit();		//Riempie un vettore fitness e ordina la popolazione


for(int j=0; j<g.Getngen(); ++j) {		//Numero di generazioni
g.Accumula();

   for(int i=0; i<1000; ++i){		

       g.Selection(r);
       g.Crossover(r);
       g.Mutate(r);
       g.newL1();
       g.Insert(i);
       i++;
      
               
}

g.Replace();
g.L1();
}

g.Print();


cout<<endl<<g.Getfit(0)<<endl;
cout<<endl;
for(int i=0; i<g.Getngenes(); ++i){

    cout<<g.Getpop(0,i)<<" ";
	
	
					}
cout<<endl<<endl;
					


return 0;
}
