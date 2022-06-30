#include "morisetti.h"


double deltay(double length, double cos){

	double y = sqrt( pow(length/2.,2) - pow(length*cos/2.,2)  );
	return y;
	
}



