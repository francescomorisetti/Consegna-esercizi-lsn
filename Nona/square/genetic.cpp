#include "genetic.h"


genetic::genetic(Random &rand) {	//Costruttore della popolazione

for(int j=0; j<npop; ++j) {
   for(int i=0; i<ngenes; ++i) {

	population[j][i]=i+1;
				}
				}
		
    random_device rd;
    mt19937 g(rd());
    				
    for (auto& row: population) { 
        std::shuffle(row.begin()+1, row.end(), g);
    				}
    

for(int i=0; i<ngenes; ++i) {

     
    posizione[i].x=rand.Rannyu(-1,1);
    posizione[i].y=rand.Rannyu(-1,1);
    
    					}
    								
}



void genetic::Fit(void) {		//Crea un vettore per il fitness e ordina la 					          generazione zero

for(int i=0; i<npop; ++i) {

    double acc=0.;
    for(int j=1; j<ngenes; ++j) {
  
        acc+=sqrt(pow(posizione[population[i][j]-1].x - posizione[population[i][j-1]-1].x ,2) +        	  pow(posizione[population[i][j]-1].y - posizione[population[i][j-1]-1].y ,2));  
  
  				}
	acc+=sqrt(pow(posizione[0].x - posizione[population[i][ngenes-1]-1].x ,2) +    	   pow(posizione[0].y - posizione[population[i][ngenes-1]-1].y ,2));  

        fitness.push_back(acc);

				}
				
for(int i=0; i<npop; ++i) {

    for(int j=i+1; j<npop; ++j) {
    
    if(fitness[j]<fitness[i]) {
    
    	swap(population[j], population[i]);
    	swap(fitness[j], fitness[i]);
    
}
}
}
}






void genetic::Selection(Random &rand) {	//Operatore di selezione

do{
    sel1=(int)(npop*pow(rand.Rannyu(),psel));
    sel2=(int)(npop*pow(rand.Rannyu(),psel));
}

while(sel1==sel2);
}




void genetic::Crossover(Random &rand) {

if(rand.Rannyu()<=pcross){

    int appo=(int)rand.Rannyu(1,ngenes);
    vector<int> appo1{};
    vector<int> appo2{};
    
    for(int i=appo; i<ngenes; ++i) {
    
    	appo1.push_back(population[sel1][i]);
    	appo2.push_back(population[sel2][i]);
    
    				       }
    
    for(int i=0; i<appo; ++i) {
    
        newgen1.push_back(population[sel1][i]);
        newgen2.push_back(population[sel2][i]);
        
    				}
    				
    for(int i=1; i<ngenes; ++i) {
        for(int j=0; j<ngenes-appo; ++j) {
            
            if(population[sel2][i]==appo1[j]) newgen1.push_back(appo1[j]); 
            if(population[sel1][i]==appo2[j]) newgen2.push_back(appo2[j]);
            
            				    }
            		          }
            		          
			}

else {
   
    for(int i=0; i<ngenes; ++i) {
        
        newgen1.push_back(population[sel1][i]);
        newgen2.push_back(population[sel2][i]);
    
                                 }
      }
}




void genetic::Mutate(Random &rand) {


if(rand.Rannyu()<=pswap) {		//Pair permutation
    
    
    int s1=(int)rand.Rannyu(1,34);
    int t1=(int)rand.Rannyu(1,34);
    int s2=(int)rand.Rannyu(1,34);
    int t2=(int)rand.Rannyu(1,34);
    		
    swap(newgen1[s1], newgen1[t1]);
    swap(newgen2[s2], newgen2[t2]);
    
    
}

		   
if(rand.Rannyu()<=pshuffle) {		//Shuffle of a sequence of length l1 and l2
  
    int s1=(int)rand.Rannyu(1,34);
    int l1=(int)rand.Rannyu(1,ngenes-s1+1);
    
    int s2=(int)rand.Rannyu(1,34);	 
    int l2=(int)rand.Rannyu(1,ngenes-s2+1);
    
    random_device rd;
    mt19937 g(rd());
    
    shuffle(newgen1.begin()+s1, newgen1.begin()+s1+l1, g);
    shuffle(newgen2.begin()+s2, newgen2.begin()+s2+l2, g);


}


if(rand.Rannyu()<=pperm) {		//Permutation of two sequences of length 3
  
    int s1=(int)rand.Rannyu(1,32); 
    int s2=(int)rand.Rannyu(1,32); 	
       	
    int p1=(int)rand.Rannyu(1,32);
    int p2=(int)rand.Rannyu(1,32);
   		
    
    swap_ranges(newgen1.begin()+s1, newgen1.begin()+s1+3, newgen1.begin()+p1);
    swap_ranges(newgen2.begin()+s2, newgen2.begin()+s2+3, newgen2.begin()+p2);
    
}


if(rand.Rannyu()<=pinvert) {		//Inversion of the order
   
    int appo1=(int)rand.Rannyu(1,34);
    int appo2=(int)rand.Rannyu(1,34);
    
    int a1=(int)rand.Rannyu(-appo1, ngenes-appo1);
    int a2=(int)rand.Rannyu(-appo2, ngenes-appo2);
    
    reverse(newgen1.begin()+appo1, newgen1.begin()+appo1+a1);
    reverse(newgen2.begin()+appo2, newgen2.begin()+appo2+a2);

}
}



void genetic::newL1(void) {		//Calcola il fitness dei due nuovi individui

     double acc1=0.;
     double acc2=0.;
     for(int i=1; i<ngenes; ++i) {

        acc1+=sqrt(pow(posizione[newgen1[i]-1].x - posizione[newgen1[i-1]-1].x ,2) +        	  pow(posizione[newgen1[i]-1].y - posizione[newgen1[i-1]-1].y ,2));  
        
        acc2+=sqrt(pow(posizione[newgen2[i]-1].x - posizione[newgen2[i-1]-1].x ,2) +        	  pow(posizione[newgen2[i]-1].y - posizione[newgen2[i-1]-1].y ,2));  
  
  				}
  				
	acc1+=sqrt(pow(posizione[0].x - posizione[newgen1[ngenes-1]-1].x ,2) +    	   pow(posizione[0].y - posizione[newgen1[ngenes-1]-1].y ,2));  

        newfit1=acc1;
        
        acc2+=sqrt(pow(posizione[0].x - posizione[newgen2[ngenes-1]-1].x ,2) +    	   pow(posizione[0].y - posizione[newgen2[ngenes-1]-1].y ,2));  

        newfit2=acc2;


}



void genetic::Insert(int index) {		//Inserisce i nuovi individui nella nuova 							  popolazione

if(Getnewfit1()<=Getfit(sel1)) { 

	newpop[index]=newgen1;
}

else newpop[index]=population[sel1];

if(Getnewfit2()<=Getfit(sel2)) { 

	newpop[index+1]=newgen2;
}

else newpop[index+1]=population[sel2];

newgen1.clear();
newgen2.clear();


}



void genetic::Replace(void) {		//Sostituisce la vecchia pop con quella nuova

    population=newpop;
    

}





void genetic::L1(void) {		//Ordina la popolazione in base al fitness

for(int i=0; i<npop; ++i) {

    double acc=0.;
    for(int j=1; j<ngenes; ++j) {
  
        acc+=sqrt(pow(posizione[population[i][j]-1].x - posizione[population[i][j-1]-1].x ,2) +        	  pow(posizione[population[i][j]-1].y - posizione[population[i][j-1]-1].y ,2));  
  
  				}
	acc+=sqrt(pow(posizione[0].x - posizione[population[i][ngenes-1]-1].x ,2) +    	   pow(posizione[0].y - posizione[population[i][ngenes-1]-1].y ,2));  

        fitness[i]=acc;

				}
				
for(int i=0; i<npop; ++i) {

    for(int j=i+1; j<npop; ++j) {
    
    if(fitness[j]<fitness[i]) {
    
    	swap(population[j], population[i]);
    	swap(fitness[j], fitness[i]);
    
}
}
}
}





double genetic::Getfit(int i) {

    return fitness[i];
    
}

double genetic::Getnewfit1(void) {

    return newfit1;
}

double genetic::Getnewfit2(void) {

    return newfit2;
}


int genetic::Getpop(int k, int h) {

    return population[k][h];

}

void genetic::Getnewpop(void) {

    for(int i=0; i<npop; ++i) {
       for(int j=0; j<ngenes; ++j) {
       
       cout<<newpop[i][j]<< " ";            
      
					}
       cout<<endl;
                               }	
}


int genetic::Getngenes(void) {

    return ngenes;

}

int genetic::Getnpop(void) {

    return npop;

}


int genetic::Getngen(void) {

    return ngen;
    
}

double genetic::Getx(int i) {

    return posizione[i].x;

}

double genetic::Gety(int i) {

    return posizione[i].y;
    
}


void genetic::Accumula(void){

double accu=0;
int n=npop/2;

for(int i=0; i<n; ++i) {

	accu+=fitness[i];
				}

stampafit.push_back(accu/n);

}



void genetic::Print(void) {

    ofstream pos, ave;
    pos.open("square.out",ios::app);
    ave.open("L.out",ios::app);
    int m;
    
    for(int i=0; i<ngenes; ++i) {
    
    m=population[0][i];
    pos<< posizione[m-1].x << " " << posizione[m-1].y<<endl;
    
    }
    pos<< posizione[0].x << " " << posizione[0].y<<endl;
    
    for(int i=0; i<ngen; ++i) {
    
    ave<<i<<" "<<stampafit[i]<<endl;
    
    }
    pos.close();
    ave.close();
    

}












