#include <algorithm>
#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include "optimize.h"

using namespace std;

int main() {

   Input();
   Initialize();
   Reset(1);

//Ricerca parametri variazionali
     
   for(int i=1; i<=SAsteps; ++i)
   {
   
   
   for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
   {
    	Reset(iblk);   //Reset block averages
   	for(int istep=1; istep <= nstep; ++istep)
    	{
    	 Move(mu, sigma);
    	 
     	 Measure(x, mu, sigma);
     	 Accumulate(); //Update block averages
    	 }
    	 
    Averages(iblk); 
   }
   glob_av=glob_av/(double)nblk;
   
   SA();
   beta=beta+0.2;
   
   }
   
   
//Simulazione con parametri fissati 


for(int iblk=1; iblk <= nblk; ++iblk) 
   {
    	Reset(iblk);   
   	for(int istep=1; istep <= nstep; ++istep)
    	{
    	 Move(mu, sigma);
    	 Fillhisto();
     	 Measure(x, mu, sigma);
     	 Accumulate(); 
    	 }
    Stampa(iblk); 
   }



return 0;

}


void Input(void)
{
  ifstream Primes, Seed, ReadInput;
  

  cout << "Variational Monte Carlo             " << endl;
  cout << "Single quantum particle in 1D" << endl;
  cout << "The program uses h=1 and m=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   Primes.open("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();
   
   
   
   Seed.open("seed.in");
   Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   Seed.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> mu;
  cout << "Mu = " << mu << endl;

  ReadInput >> sigma;
  cout << "Sigma = " << sigma << endl;

  ReadInput >> x;
  cout << "Initial position = " << x << endl;

  ReadInput >> delta;
  cout << "Sampling distance  = " << delta << endl << endl;
  
  ReadInput >> beta;
  
  ReadInput >> SAsteps;
  
  ReadInput >> eps;
    
  ReadInput >> nblk;

  ReadInput >> nstep;
  
  ReadInput >> nbins;

  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

}


void Initialize(void)
{


    	Reset(1);   //Reset block averages
   	for(int istep=1; istep <= nstep; ++istep)
    	{
    	 Move(mu, sigma);
     	 Measure(x, mu, sigma);
     	 Accumulate(); //Update block averages
    	 }
    Averages(1);   //Print results for current block
    Print();


}


void Move(double m, double s)
{
  
  double A, x_old, x_new;
  
  x_old = x;
  x_new = x + rnd.Rannyu(-delta,delta);
  A = min(1., pow(Psi(x_new, m, s),2)/pow(Psi(x_old, m, s),2));
  
  if(rnd.Rannyu()<=A) 
  {
  x=x_new;
  accepted++;
  }  
  
  attempted++;	
}

void SA(void)
{

   double P, mu_new, sigma_new;
   
   
   mu_new = mu + rnd.Rannyu(-eps,eps);
   sigma_new = sigma + rnd.Rannyu(-eps,eps);
   
  
   for(int iblk=1; iblk <= nblk; ++iblk) 
   {
    	
    	if(iblk == 1)
   {
       glob_new = 0;
       glob_new2 = 0;
   }

   
   blk_av = 0;
   
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
    	
   	for(int istep=1; istep <= nstep; ++istep)
    	{
    	 Move(mu_new, sigma_new);
     	 Measure(x, mu_new, sigma_new);
     	 Accumulate(); 
    	 }
    stima_H = blk_av/blk_norm; 
    glob_new  += stima_H;
    glob_new2 += stima_H*stima_H;
    err_new=Error(glob_new,glob_new2, iblk);    
    
   }
   glob_new=glob_new/(double)nblk;
   
   
   P = exp(-beta*(glob_new - glob_av));
   
   if(glob_new<=glob_av) 
   {
   mu = mu_new;
   sigma = sigma_new;
   glob_av=glob_new;
   err_H=err_new;
   
   }
   
   else if(rnd.Rannyu()<=P) 
   {
   mu = mu_new;
   sigma = sigma_new;
   glob_av=glob_new;
   err_H=err_new;
   }
   
   Print();
   
}

void Fillhisto(void)
{
    
    ofstream funz;
    funz.open("density.out",ios::app);
    
    funz << x << " " << pow(Psi(x, mu, sigma),2) <<endl;
    
    funz.close();
       	
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       glob_av = 0;
       glob_av2 = 0;
   }

   
   
   blk_av = 0;
   
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Measure(double x, double m, double s) 
{
    E=Ene(x, m, s);

}

void Accumulate(void) //Update block averages
{

   
     blk_av = blk_av + E;
   
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    
    stima_H = blk_av/blk_norm; 
    glob_av  += stima_H;
    glob_av2 += stima_H*stima_H;
    err_H=Error(glob_av,glob_av2,iblk);    
    

     cout << "----------------------------" << endl << endl;
}

void Stampa(int iblk)
{
   ofstream energia;
   energia.open("ene.out",ios::app);

   stima_H = blk_av/blk_norm; 
   glob_av  += stima_H;
   glob_av2 += stima_H*stima_H;
   err_H=Error(glob_av,glob_av2,iblk);    

   energia << iblk << " " << glob_av/(double)iblk << " " << err_H <<endl;
   energia.close();

}

void Print(void) 		//Stampa i risultati per ogni temperatura
{

    ofstream out;
    
    out.open("Variational.out",ios::app);
    out << beta <<  " " << glob_av << " " << err_H << " "<< mu << " " << sigma << endl;
    out.close();

}


double Psi(double x, double m, double s) {


    return exp( -pow(x-m,2)/(2.*pow(s,2)) ) + exp( -pow(x+m,2)/(2.*pow(s,2)) ); 


}

double Ene(double x, double m, double s) {

    return  -0.5*(exp( -pow(x-m,2)/(2.*pow(s,2)) )*pow(m-x,2)/pow(s,4) + exp( -pow(x+m,2)/(2.*pow(s,2)) )*pow(m+x,2)/pow(s,4) - Psi(x, m, s)/pow(s,2))/Psi(x, m, s) + pow(x,4)-2.5*pow(x,2);

}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}


