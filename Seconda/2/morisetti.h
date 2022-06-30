#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <assert.h>
#include "random.h"

using namespace std;

struct coordinate {

	double x;
	double y;
	double z;
		};
		
double moduloquadro(double, double, double);


#ifndef __funzioni_h__
#define __funzioni_h__

class FunzioneBase{

public:

  virtual double Eval(double x) const = 0;

};

class Coseno: public FunzioneBase{

public:

  Coseno() {
	  m_a=0;
  };

  Coseno(double a) {
	  m_a=a;
  };

  ~Coseno();
  virtual double Eval(double x) const{return (cos(m_a*x));};
  void SetA(double a) {m_a = a;};
  double GetA() const {return m_a;};
  

private:

  double m_a;

};

#endif

#ifndef __randomwalk_h__
#define __randomwalk_h__


class randomwalk {

 public:

  randomwalk(unsigned int N, Random &r);
  //continuo(unsigned int N, Random &r); 

  ~randomwalk();
  //~continuo();

	
  int Getx(unsigned int) ;
  int Gety(unsigned int) ;
  int Getz(unsigned int) ;
  
  coordinate& operator[] (unsigned int);

 private:

	unsigned int length;
	coordinate* m_v;

};

#endif // __randomwalk_h__


#ifndef __continuo_h__
#define __continuo_h__


class continuo {

 public:

  continuo(unsigned int N, Random &r);
  

  ~continuo();
  

	
  double Getx(unsigned int) ;
  double Gety(unsigned int) ;
  double Getz(unsigned int) ;
  
  coordinate& operator[] (unsigned int);

 private:

	unsigned int length;
	coordinate* m_c;

};

#endif // __continuo_h__


