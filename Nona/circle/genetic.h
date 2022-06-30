#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <iterator>
#include <random>
#include "random.h"

using namespace std;

struct coordinate {

	double x;
	double y;
		};
		
	

#ifndef __genetic_h__
#define __genetic_h__




class genetic {

 public:

 genetic(Random &rand); 	//Costruisce la popolazione

 
 void Fit(void);
 void L1(void);
 void newL1(void);
 void Insert(int);
 void Selection(Random &rand);
 void Mutate(Random &rand);
 void Crossover(Random &rand);
 void Replace(void);
 
 int Getpop(int, int);
 int Getnpop(void);
 int Getngen(void);
 void Getnewpop(void);
 int Getngenes(void);	
 double Getfit(int);
 double Getnewfit1(void);
 double Getnewfit2(void);
 
 double Getx(int);
 double Gety(int);
 
 void Print(void);
 void Accumula(void);
 

 private:

 int ngenes = 34;
 int npop = 1000;
 int ngen = 500;
  
 int sel1, sel2;
 double pcross=0.6;
 int psel=5;
 double pswap=0.2;
 double pshuffle=0.2;
 double pperm=0.1;
 double pinvert=0.1;
 
 vector<int> newgen1{};
 vector<int> newgen2{};
 
 vector<vector<int>> population{npop, vector<int>(ngenes)};
 vector<vector<int>> newpop{npop, vector<int>(ngenes)};
 vector<coordinate> posizione{ngenes};
 vector<double> fitness{};
 vector<double> stampafit{};
 double newfit1;
 double newfit2;
 
};

#endif // __genetic_h__






