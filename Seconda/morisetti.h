#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;


double deltay(double, double);

double traj(int, int); 


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
