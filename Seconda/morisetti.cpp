#include "morisetti.h"


double deltay(double length, double cos){

	double y = sqrt( pow(length/2.,2) - pow(length*cos/2.,2)  );
	return y;
	
}

double traj(int length, int dim) {

double* pos=new double[length][dim];
double step;


for(int i=0; i<length; i++) {
	
	pos[i][:]=0;
	step=rnd.Rannyu();
	
	if (step<(1/3)) {
		if (step%2==0) pos[i][0]+=1;
		else if (step%2!=0) pos[i][0]+=-1;
	}
	
	if (step>=(1/3) and step<(2/3)) {
		if (step%2==0) pos[i][1]+=1;
		else if (step%2!=0) pos[i][1]+=-1;	
	  
	}
	
	if (step>=(2/3)) {
		if (step%2==0) pos[i][2]+=1;
		else if (step%2!=0) pos[i][2]+=-1;

	}
				

return pos;




}

