#include "morisetti.h"
//#define NDEBUG 


continuo::continuo(unsigned int N, Random &r) {	//costruttore rw continuo

	length = N;
	m_c = new coordinate[length];
	double theta;
	double phi;
	
	theta=r.Theta();
  	phi=r.Phi();
  	
	m_c[0].x=sin(theta)*cos(phi);
  	m_c[0].y=sin(theta)*sin(phi);
  	m_c[0].z=cos(theta);
			
  	
  	for ( int i = 1 ; i < length ; i++) {
  	
  		theta=r.Theta();
  		phi=r.Phi();
  		
  		m_c[i].x=m_c[i-1].x + sin(theta)*cos(phi);
  		m_c[i].y=m_c[i-1].y + sin(theta)*sin(phi);
  		m_c[i].z=m_c[i-1].z + cos(theta);
  		
  		
  		
  						}

}


continuo::~continuo() {			//distruttore rw discreto
  
  delete[]m_c;
}


double continuo::Getx(unsigned int i)  {

  if ( i<length ) {
    return m_c[i].x;
  } else {
    cerr << "Errore: indice " << i << ", dimensione " << length << endl;
    exit (-1);

  }
}

double continuo::Gety(unsigned int i)  {

  if ( i<length ) {
    return m_c[i].y;
  } else {
    cerr << "Errore: indice " << i << ", dimensione " << length << endl;
    exit (-1);

  }
}

double continuo::Getz(unsigned int i)  {

  if ( i<length ) {
    return m_c[i].z;
  } else {
    cerr << "Errore: indice " << i << ", dimensione " << length << endl;
    exit (-1);

  }
}

coordinate& continuo::operator[] (unsigned int i) {

  assert( ( length > i ) && "Errore : l'indice e' troppo grande");
  if ( i<length ) {
    return m_c[i];
} else {
    cerr << "Errore: indice " << i << ", dimensione " << length << endl;
    exit (-1);
  }
}



randomwalk::randomwalk(unsigned int N, Random &r) {	//costruttore rw discreto
   
  length = N;
  m_v = new coordinate[length];
  double pos;
  
  	
  for ( int i = 0 ; i < length ; i++) {
  	m_v[i].x=0;
  	m_v[i].y=0;
  	m_v[i].z=0;
  					}
  	pos=r.Rannyu(0,3);
  	
	if (pos<1) {
		if (r.Rannyu()<0.5) m_v[0].x=1;
		else m_v[0].x=-1;
	}
	
	if (pos>=1 and pos<2) {
		if (r.Rannyu()<0.5) m_v[0].y=1;
		else m_v[0].y=-1;	
	  
	}
	
	if (pos>=2) {
		if (r.Rannyu()<0.5) m_v[0].z=1;
		else m_v[0].z=-1;

	}
  
  
  for ( int i = 1 ; i < length ; i++) {
  
	pos=r.Rannyu(0,3);	
	
	if (pos<1) {
	
		m_v[i].y=m_v[i-1].y;
		m_v[i].z=m_v[i-1].z;
		
		if (r.Rannyu()<0.5) m_v[i].x=m_v[i-1].x + 1;
		else m_v[i].x=m_v[i-1].x - 1;
	}
	
	if (pos>=1 and pos<2) {
		
	  	m_v[i].x=m_v[i-1].x;
		m_v[i].z=m_v[i-1].z;
		
		if (r.Rannyu()<0.5) m_v[i].y=m_v[i-1].y + 1;
		else m_v[i].y=m_v[i-1].y - 1;
	}
	
	if (pos>=2) {
	
		m_v[i].y=m_v[i-1].y;
		m_v[i].x=m_v[i-1].x;
		
		if (r.Rannyu()<0.5) m_v[i].z=m_v[i-1].z + 1;
		else m_v[i].z=m_v[i-1].z - 1;

	}
	}			
}

randomwalk::~randomwalk() {			//distruttore
  
  delete[]m_v;
}

int randomwalk::Getx(unsigned int i)  {

  if ( i<length ) {
    return m_v[i].x;
  } else {
    cerr << "Errore: indice " << i << ", dimensione " << length << endl;
    exit (-1);

  }
}

int randomwalk::Gety(unsigned int i)  {

  if ( i<length ) {
    return m_v[i].y;
  } else {
    cerr << "Errore: indice " << i << ", dimensione " << length << endl;
    exit (-1);

  }
}

int randomwalk::Getz(unsigned int i)  {

  if ( i<length ) {
    return m_v[i].z;
  } else {
    cerr << "Errore: indice " << i << ", dimensione " << length << endl;
    exit (-1);

  }
}

coordinate& randomwalk::operator[] (unsigned int i) {

  assert( ( length > i ) && "Errore : l'indice e' troppo grande");
  if ( i<length ) {
    return m_v[i];
} else {
    cerr << "Errore: indice " << i << ", dimensione " << length << endl;
    exit (-1);
  }
}

double moduloquadro(double x, double y, double z){

return pow(x,2)+pow(y,2)+pow(z,2);


}
