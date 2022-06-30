

#ifndef __ISING__
#define __ISING__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables

double mu, sigma;
double x;
double delta;
double beta;
double SAsteps;
double eps;



// averages
double blk_av,blk_norm,accepted,attempted;
double glob_av,glob_av2, glob_new, glob_new2;
double stima_H;
double err_H, err_new;


// simulation
int nstep, nblk, nbins;
double E;

//functions
void Input(void);
void Initialize(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Fillhisto(void);
void Print(void);
void Stampa(int);

void Move(double, double);
void SA(void);

void Measure(double, double, double);
double Psi(double, double, double);
double Ene(double, double, double);
double Error(double,double,int);

#endif


