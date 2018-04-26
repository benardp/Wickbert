/**
 * Implementation of a specular implicit 
 * @file Icosahedral.cpp
 * @date 16 March. 2006
 * @author Matei N. Stroila
 * @remarks
 */

#include "Icosahedral.h"

static double T = 1.618; //golden ratio 

REGISTER_IMPLICIT(Icosahedral,"Icosahedral");

Icosahedral::Icosahedral(void)
{
}

double Icosahedral::proc(const gmVector3 & x)
{
	return 2 - (cos(x[0] + T*x[1]) + cos(x[0] - T*x[1]) + cos(x[1] + T*x[2]) + cos(x[1] - T*x[2]) + cos(x[2] - T*x[0]) + cos(x[2] + T*x[0]));
}


gmVector3 Icosahedral::grad(const gmVector3 & x)
{
	gmVector3 gradient;
	gradient[0] = sin(x[0] + T*x[1]) + sin(x[0] - T*x[1]) - T * sin(x[2] - T*x[0]) + T * sin(x[2] + T*x[0]);
	gradient[1] = T * sin(x[0] + T*x[1]) - T * sin(x[0] - T*x[1]) + sin(x[1] + T*x[2]) + sin(x[1] - T*x[2]);
	gradient[2] = T * sin(x[1] + T*x[2]) - T * sin(x[1] - T*x[2]) + sin(x[2] - T*x[0]) + sin(x[2] + T*x[0]);
	return gradient;
}


gmMatrix3 Icosahedral::hess(const gmVector3 & x) 
{
	double hxx = cos(x[0] + T*x[1]) + cos(x[0] - T*x[1]) + T * T * cos(x[2] - T*x[0]) + T * T * cos(x[2] + T*x[0]);
	double hxy = T * cos(x[0] + T*x[1]) - T * cos(x[0] - T*x[1]);
	double hxz = - T * cos(x[2] - T*x[0]) + T * cos(x[2] + T*x[0]);
	double hyy =  T * T * cos(x[0] + T*x[1]) + T * T * cos(x[0] - T*x[1]) + cos(x[1] + T*x[2]) + cos(x[1] - T*x[2]);
	double hyz =  T * cos(x[1] + T*x[2]) - T * cos(x[1] - T*x[2]);
	double hzz = T * T * cos(x[1] + T*x[2]) + T * T * cos(x[1] - T*x[2]) + cos(x[2] - T*x[0]) + cos(x[2] + T*x[0]);
	gmMatrix3 hessian(hxx,  hxy,  hxz,
                    hxy,  hyy,  hyz,
                    hxz,  hyz,  hzz);
	return hessian;
}
