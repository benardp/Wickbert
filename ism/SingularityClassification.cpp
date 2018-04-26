/**
 * Implementation of a class used to classify the critical points of a feature contour.
 * @file SingularityClassification.h
 * @date 10/14/2005
 * @author Matei N. Stroila
 * @remarks 
 */

#include "Surface/Implicit/Implicit.h"
#include "Newton.h"
#include "SingularityClassification.h"

SingularityClassification::SingularityClassification(Implicit *myF)
{
	F = myF;
}

SingularityType SingularityClassification::getSingularityType(gmVector3& myP)
{
	p = myP;
	computeSigma();
	if(gmIsZero(Sigma)) return DEGENERATE;
	if(Sigma > 0) return EXTREMUM;
	if(Sigma < 0) return SADDLE;
	return DEGENERATE;
}

void SingularityClassification::computeSigma()
{
	double Fx, Fy, Fxx, Fyy, Fyz, Fxy, Fxz, Fxxz, Fxzz, Fyzz, Fyyz, Fxyz, Fzzz;
	Fx = F->Fx(p);
	Fy = F->Fx(p);
	Fxx = F->Fxx(p);
	Fxy = F->Fxy(p);
	Fxz = F->Fxz(p);
	Fyy = F->Fyy(p);
	Fyz = F->Fyz(p);
	Fxxz = F->Fxxi(p,2);
	Fxzz = F->Fxzi(p,2);
	Fyzz = F->Fyzi(p,2);
	Fyyz = F->Fyyi(p,2);
	Fxyz = F->Fxyi(p,2);
	Fzzz = F->Fzzi(p,2);
	if(!gmIsZero(Fx)){
		Sigma = - pow(Fx,2) * pow(Fxzz,2) * pow(Fy,2) + 2 * pow(Fx,3) * Fxzz * Fy * Fyzz
			- pow(Fx,4) * pow(Fyzz,2) - 2 * pow(Fx,3) *  Fxyz * Fy * Fzzz
			+ 2 * pow(Fx,2) *  Fxy * Fxz * Fy * Fzzz +  pow(Fx,2) * Fxxz * pow(Fy,2) * Fzzz
			- Fx * Fxx * Fxz * pow(Fy,2) * Fzzz -  pow(Fx,3) * Fxz * Fyy * Fzzz
			+ pow(Fx,4) * Fyyz * Fzzz;
		return;
	}
	if(!gmIsZero(Fy)){
		Sigma = - pow(Fxzz,2) * pow(Fy,4) + 2 * Fx * Fxzz * pow(Fy,3) * Fyzz
		- pow(Fx,2) * pow(Fy,2) * pow(Fyzz,2) - 2 * Fx *  Fxyz * pow(Fy,3) * Fzzz 
		+ Fxxz *  pow(Fy,4) * Fzzz +  pow(Fx,2) * pow(Fy,2) * Fyyz * Fzzz
		+ 2 * Fx * Fxy * pow(Fy,2) * Fyz * Fzzz -  Fxx * pow(Fy,3) * Fyz * Fzzz
		- pow(Fx,2) * Fy * Fyy * Fyz * Fzzz;
	}
	
}
