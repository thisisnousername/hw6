/*
author: markus
  date: 2015-12-03
*/

#include <iostream>
#include <cmath>

using namespace std;

void f(double* y, double* k, double* k1, double* k2, double* k3, double* h, const double a, const double b, const double c, const double dt){

	k[0] = a*((y[1]+dt*h[0]*k1[1]+dt*h[1]*k2[1]+dt*h[2]*k3[1])-(y[0]+dt*h[0]*k1[0]+dt*h[1]*k2[0]+dt*h[2]*k3[0]));
	k[1] = (y[0]+dt*h[0]*k1[0]+dt*h[1]*k2[0]+dt*h[2]*k3[0])*(b-(y[2]+dt*h[0]*k1[2]+dt*h[1]*k2[2]+dt*h[2]*k3[2]))+(y[1]+dt*h[0]*k1[1]+dt*h[1]*k2[1]+dt*h[2]*k3[1]);
	k[2] = (y[0]+dt*h[0]*k1[0]+dt*h[1]*k2[0]+dt*h[2]*k3[0])*(y[1]+dt*h[0]*k1[1]+dt*h[1]*k2[1]+dt*h[2]*k3[1])-c*(y[2]+dt*h[0]*k1[2]+dt*h[1]*k2[2]+dt*h[2]*k3[2]);

}

int main(){

	const int N=10000;
	const double dt=100.0/N*1.0;

	const double a=10;
	const double b=28;
	const double c=8/3;

	double y[3];
	double k1[3];
	double k2[3];
	double k3[3];
	double k4[3];
	double h[3];

	y[0] = 1;
	y[1] = 1;
	y[2] = 1;

	cout << 0 << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] <<endl;

	for(int i=0; i<N; i++){

		h[0]=0; h[1]=0; h[3]=0;
		f(y, k1, k1, k2, k3, h, a, b, c, dt);

		h[0]=0.5; h[1]=0; h[3]=0;
		f(y, k2, k1, k2, k3, h, a, b, c, dt);

		h[0]=0; h[1]=0.5; h[3]=0;
		f(y, k3, k1, k2, k3, h, a, b, c, dt);

		h[0]=0; h[1]=0; h[3]=1;
		f(y, k4, k1, k2, k3, h, a, b, c, dt);

		y[0] = y[0] +dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
		y[1] = y[1] +dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
		y[2] = y[2] +dt/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);


	cout << i*dt << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] <<endl;
	}

	return 0;

}
