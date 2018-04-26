/*
 *  RotatedImplicit.cpp
 *  gps
 *
 *  Created by Matei Stroila on 10/13/05.
 *  
 */

#include "RotatedImplicit.h"

REGISTER_IMPLICIT(RotatedImplicit,"RotatedImplicit");

RotatedImplicit::RotatedImplicit(Implicit *F, gmMatrix3& R)
{
	myF = F;
	myR = R;
	rotHasChanged = false;
}

RotatedImplicit::RotatedImplicit(void)
{
	myF = NULL;
	myR = myR.identity();
	new SurfImpRefParam(this,&myF,"<empty>","implicit","Implicit",
		"Implicit whom silhouette we want");
	new SurfParamBool(this,&rotHasChanged, true,"Rotation Has Changed","Rotation Has Changed","Set True if want to update the Rotation");
}

void RotatedImplicit::setF(Implicit *F)
{
	myF = F;
	rotHasChanged = false;
}

void RotatedImplicit::setRotation(gmMatrix3& R)
{
	myR = R;
	rotHasChanged = false;
}

double RotatedImplicit::proc(const gmVector3 & x)
{
	if(rotHasChanged)
	{	
		checkRotation();
	}
	if(!myF) return 0; 
	return myF->proc(myR * x);
	
}

Intervald RotatedImplicit::proc(const Box<double>&  b)
{
	if(rotHasChanged)
	{	
		checkRotation();
	}
	if(!myF) return Intervald(0);
	Box<double> X(3);
	for(unsigned int i = 0; i<3; i++)
		for(unsigned int j = 0; j<3; j++) 
			X[i] += myR[i][j] * b[j];
	
	return myF->proc(X);
	
}

gmVector3 RotatedImplicit::grad(const gmVector3 & x)
{
	if(rotHasChanged)
	{	
		checkRotation();
	}
	if(!myF) return gmVector3(-0.1,-0.1,-0.1); 
	return myR.transpose() * (myF->grad( myR * x));
}



gmMatrix3 RotatedImplicit::hess(const gmVector3 & x)
{
	if(rotHasChanged)
	{	
		checkRotation();
	}
	if(!myF) return myR; 
	return myR.transpose() * (myF->hess( myR * x)) * myR;
}

void RotatedImplicit::checkRotation(void) 
{	
	GLdouble modelview[16];
	glGetDoublev(GL_MODELVIEW_MATRIX,modelview);
	gmMatrix4 _m(modelview[0], modelview[4], modelview[8], modelview[12],
		modelview[1], modelview[5], modelview[9], modelview[13],
		modelview[2], modelview[6], modelview[10], modelview[14],
		modelview[3], modelview[7], modelview[11], modelview[15]);
	for(unsigned int i = 0; i < 3; i++)
		for(unsigned int j = 0; j < 3; j++)
		myR[i][j] = _m[i][j];
	
	rotHasChanged = false;
}



 
 