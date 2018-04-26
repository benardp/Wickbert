/*
 *  ZeroGaussianCurvature.h
 *  gps
 *
 *  Created by Matei Stroila on 9/27/05.
 *  
 */

#ifndef ZeroGaussianCurvature_h
#define ZeroGaussianCurvature_h

#include "Surface/Implicit/Implicit.h"
#include "libgm/gm.h"

class ZeroGaussianCurvature : public Implicit
{
public:
	ZeroGaussianCurvature();
	ZeroGaussianCurvature(Implicit *F); ///< Explicit constructor.
	
	double proc(const gmVector3 & x);
	Intervald proc(const Box<double>&  b);
	//gmVector3 grad(const gmVector3 & x);
	
	~ZeroGaussianCurvature(){}; 
	MAKE_NAME();

 private:
		
	Implicit *myF; ///< The input implicit.
	
};

#endif

