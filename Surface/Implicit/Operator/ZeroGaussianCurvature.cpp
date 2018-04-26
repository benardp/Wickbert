/*
 *  ZeroGaussianCurvature.cpp
 *  gps
 *
 *  Created by Matei Stroila on 9/27/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "ZeroGaussianCurvature.h"

REGISTER_IMPLICIT(ZeroGaussianCurvature,"ZeroGaussianCurvature");

ZeroGaussianCurvature::ZeroGaussianCurvature(Implicit *F)
{
	myF = F;
}

ZeroGaussianCurvature::ZeroGaussianCurvature(void)
{
	myF = NULL;
	new SurfImpRefParam(this,&myF,"<empty>","implicit","Implicit",
						"Implicit whom silhouette we want");
}

double ZeroGaussianCurvature::proc(const gmVector3 & x)
{
	return myF->numeratorGaussianCurvature(x);
}

//TO DO: implement this
Intervald ZeroGaussianCurvature::proc(const Box<double>&  b)
{
	return Intervald();
}

//gmVector3 ZeroGaussianCurvature::grad(const gmVector3 & x)
//{
	//return myF->gradGaussianCurvature(x);
//}

 
 