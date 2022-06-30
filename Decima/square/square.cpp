#include "genetic.h"



int main() {


	MPI_Init(&argc,&argv);
	int node_rank, nodes_number;

	MPI_Comm_rank(MPI_COMM_WORLD, &node_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nodes_number);

	Random r;
	
	int seed[4];
   int* p1 = new int[nodes_number]();
   int* p2 = new int[nodes_number]();
   ifstream Primes("Primes");
   
   for(int irank=0; irank<nodes_number; irank++){
   
   	if (Primes.is_open()){
           Primes >> p1[irank] >> p2[irank] ;
      } else cerr << "PROBLEM: Unable to open Primes" << endl;
      }
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            r.SetRandom(seed,p1[node_rank],p2[node_rank]);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	
	Random r;
	r=r.inizializza(r);

	genetic g(r);		//Crea la popolazione iniziale e le posizioni spaziali
	g.Fit();		//Riempie un vettore fitness e ordina la popolazione

vector<double> x(50);


for(int j=1; j<=g.Getngen(); ++j) {		//Numero di generazioni
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
	
	if(j%50==0){
	
		for(int node=0; node<nodes_number; ++node){
	           
	           for(int i=0; i<g.Getngenes(); ++i) {	//Vettore da inviare
			x[i]=g.Getpop[0][i];
							}
		MPI_Bcast(&x.front(), x.size(), MPI_DOUBLE, nodo, MPI_COMM_WORLD);	
			
			 				}
			}
	g.L1();
			
			
}

	g.Print(node_rank);
	
	
	
					
MPI_Finalize();

return 0;
}











