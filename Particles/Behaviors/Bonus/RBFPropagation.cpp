//RBF Propagation
//Mike Mullan, John C. Hart

#include "RBFPropagation.h"
#include "ParticlePosition.h"

#include <iostream>
using namespace std;

REGISTER_PARTICLESTUFF(RBFPropagation,"Behavior:RBFPropagation");

using namespace TNT;

void RBFPropagation::attachAttributes()
{
	attachAttribute(imp_int,std::string("ImplicitInterrogator"));
	attachAttribute(p_orient,std::string("ParticleOrientation"));
}


RBFPropagation::RBFPropagation(Particles *ps, const std::string& name)
: ParticleBehavior(ps, std::string("RBFPropagation"))
{
	speedfac = 0.0;
	//constfac = 1.0;
	constfac = 0.0;
	gaussfac = 0.0;
	meanfac = 0.0;
	attachAttributes();
}


double RBFPropagation::speed(gmVector3 point)
{
	double s;

	s = constfac;

	if (!imp_int)
		return s;
	RBFInterpolant *rbf = dynamic_cast<RBFInterpolant *>(imp_int->getImplicit());
	if (!rbf)
		return s;
	
	//Speed based on the distance to a target surface
	if(rbf->morph)
	{
		s += rbf->targetShape(point) ;
		
	
	}
	//Speed based on curvature
	else
	{
		s += gaussfac*rbf->gaussianCurvature(point);
		s += meanfac*rbf->meanCurvature(point);
	}
	
	s *= speedfac;
	
	return s;
	
}


void RBFPropagation::applyConstraint()
{
	if (!imp_int) return;

	// imp_int->getImplicit()->interpolate(ps,flexible);
	
	RBFInterpolant *rbf = dynamic_cast<RBFInterpolant *>(imp_int->getImplicit());

	if (!rbf) return;

	if (speedfac == 0.0)
		return;
	
	 //Move the rbf dipoles in their normal direction by the speed function, and then adjust the normals to reflect the new shape of the rbf.
	
	double eps = -rbf->centers[0].h;
	rbf->moving = 1;	
		
	double speedTol = .01;
	int numOnSurface = 0;

	//Move by the speed function in the direction that reduces the speed function (for morphing) 
	for(unsigned int i = 0; i < rbf->normals.size(); i++)
	{
		double s = speed(  rbf->centers[i].c );
		rbf->centers[i].c += s * rbf->normals[ i ];
		double plus = speed(rbf->centers[i].c);
		rbf->centers[i].c -= 2 * s * rbf->normals[ i ];
		double minus = speed(rbf->centers[i].c);
		
		if(fabs(plus) < fabs(s) )
			rbf->centers[i].c += 2 * s * rbf->normals[ i ];
		
		if(fabs(s) < fabs(speedTol) )
			numOnSurface++;
	}

	
	if(numOnSurface == rbf->normals.size() )
	{
		printf("Tuning\n");
		tune(rbf);
	}
	
//	adaptivity(rbf);

	//for(unsigned int i = 0; i < rbf->normals.size(); i++)
	//		rbf->centers[i].c += speed(  rbf->centers[i].c ) * rbf->normals[ i ];

	

	//The normal vectors of RBF dipoles will rotate with the propagation of the surface
	if(rbf->dipoles)
		propagateNormals(rbf);
	else
	{
		for(unsigned int i = 0; i < rbf->normals.size(); i++)
		{
			rbf->normals[i] = rbf->grad(rbf->centers[i].c);
			rbf->normals[i].normalize();
		}	
	}

	rbf->updateRBF();
	
}

void RBFPropagation::integrate()
{
}

//Rotates the dipole normals
void RBFPropagation::propagateNormals(RBFInterpolant* rbf)
{
	int numDipoles = rbf->normals.size() / 2;
    double eps = -rbf->centers[0].h;
	for(int i = 0; i < numDipoles; i++)
	{
		gmVector3 point = rbf->centers[2 * i].c + eps * rbf->normals[2 * i];
		
		gmMatrix3 T;
		T = T.identity();
		gmVector3 N; 
		N = rbf->normals[2 * i];
		
		for(int row = 0; row < 3; row++)
			for(int col = 0; col < 3; col++)
			   	T[row][col] = T[row][col] - (N[row] * N[col]);
				
		//gmMatrix3 hessian = rbf->computeHessian(point);
		gmMatrix3 hessian = rbf->hess(point);
		
		gmVector3 nDot = speed(point) * (T * hessian * N) + (rbf->grad(point)).length() * (T * computeGradF(point)); //eq 8
		

		//	gmVector3 nDot =  -(rbf->grad(point)).length() * (T * computeGradF(point));
		rbf->normals[2 * i] += nDot;
		rbf->normals[2 * i + 1] += nDot;
		rbf->normals[2 * i].normalize();
		rbf->normals[2 * i + 1].normalize();
		rbf->centers[i * 2 + 1].c = rbf->centers[i * 2].c + 2 * eps * rbf->normals[2 * i];

		
	}
}

//The gradient of the speed function
gmVector3 RBFPropagation::computeGradF(gmVector3 point)
{
	double eps = .0001;
	gmVector3 gradF;
	gradF[0] = ( speed( gmVector3(point[0] + eps, point[1], point[2] )) - speed( gmVector3(point[0] - eps, point[1], point[2] )) ) / ( 2 * eps );
	gradF[1] = ( speed( gmVector3(point[0], point[1] + eps, point[2] )) - speed( gmVector3(point[0], point[1] - eps, point[2] )) ) / ( 2 * eps );
	gradF[2] = ( speed( gmVector3(point[0], point[1], point[2] + eps )) - speed( gmVector3(point[0], point[1], point[2] - eps )) ) / ( 2 * eps );
	return gradF;
}


double RBFPropagation::targetShape(gmVector3 point)
{
	//return - ( pow(point[0], 4) -  pow(point[0], 2) +  pow(point[1], 4) -  pow(point[1], 2) +  pow(point[2], 4) -  pow(point[2], 2) - .5);
	
	return  - ( point[0] * point[0] + point[1] * point[1] + point[2] * point[2] - 1 );
}

void RBFPropagation::adaptivity(RBFInterpolant* rbf)
{
	for(unsigned int i = 0; i < rbf->normals.size(); i++)
	{
		if( fabs( rbf->meanCurvature(rbf->centers[i].c) ) > 100)
		{
			rbf->centers.erase( rbf->centers.begin() + i);
			rbf->normals.erase( rbf->normals.begin() + i);
		}	
	}

	ParticlePosition* positions;
	positions = ps->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
	assert(positions);

	double distance;
	double minDistance = 1000;
	double cutoff = 1.0;
	gmVector3 newNormal;

	for(unsigned int i = 0; i < ps->size(); i++)
	{
		for(unsigned int j = 0; j < rbf->normals.size(); j++)
		{
			gmVector3 centerToParticle = rbf->centers[j].c - positions->x[i];
			distance = centerToParticle.length();
			
			if(distance < minDistance)
				minDistance = distance;
		}
	
		if(minDistance > cutoff)
		{
			if(rbf->dipoles)
			{
				
				newNormal = rbf->grad(positions->x[i]);
				newNormal.normalize();
				rbf->normals.push_back(newNormal);
				rbf->normals.push_back(newNormal);
				rbf->centers.push_back(RBFPoint(positions->x[i] - epsilon * newNormal, -epsilon));
				rbf->centers.push_back(RBFPoint(positions->x[i] + epsilon * newNormal, epsilon));
			
				
			
			}
			else
			{
				RBFPoint temp = rbf->centers[rbf->normals.size()];
				rbf->centers[rbf->normals.size()] = RBFPoint(positions->x[i], 0);
                rbf->centers.push_back(temp);
				rbf->updateRBF();
				newNormal = rbf->grad(positions->x[i]);
				newNormal.normalize();
				rbf->normals.push_back(newNormal);
			}
				
			break;
		}
		
		minDistance = 1000;

	}

}

void RBFPropagation::tune(RBFInterpolant* rbf)
{
	double speedTol = .25;

	ParticlePosition* positions;
	positions = ps->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
	assert(positions);

	gmVector3 newNormal;

	for(unsigned int i = 0; i < ps->size(); i++)
	{
		if(rbf->targetMesh.proc(positions->x[i]) > speedTol)
		{
			if(rbf->dipoles)
			{
				
				newNormal = rbf->grad(positions->x[i]);
				newNormal.normalize();
				rbf->normals.push_back(newNormal);
				rbf->normals.push_back(newNormal);
				rbf->centers.push_back(RBFPoint(positions->x[i] - epsilon * newNormal, -epsilon));
				rbf->centers.push_back(RBFPoint(positions->x[i] + epsilon * newNormal, epsilon));
			
				
			
			}
			else
			{
				RBFPoint temp = rbf->centers[rbf->normals.size()];
				rbf->centers[rbf->normals.size()] = RBFPoint(positions->x[i], 0);
                rbf->centers.push_back(temp);
				rbf->updateRBF();
				newNormal = rbf->grad(positions->x[i]);
				newNormal.normalize();
				rbf->normals.push_back(newNormal);
			}
		
		}
	}
}
