/**
* @file SPHTension.cpp
* @author Yuan Zhou
*/

#include "SPHTension.h"

#define mass_default 1.0 
#define kernel_h_default 1.0 // 0.11//1.0
#define tensionC_default 1.0 // 0.1   // tension coefficient
#define thresh_default 5 //1.0  // 1e-4   // normal threshold for judge surface particles

REGISTER_PARTICLESTUFF(SPHTension,"Behavior:SPHTension");

SPHTension::SPHTension(Particles *ps) 
: ParticleBehavior(ps, std::string("SPHTension") )
{
	density = NULL;
	acceleration = NULL;
	plocality=NULL;
	mass = mass_default;
	kernel_h = kernel_h_default;
	tensionC = tensionC_default;
	normal_thresh = thresh_default;

	h2 = kernel_h*kernel_h;
	h9 = h2*h2*h2*h2*kernel_h;
	coeff = 315.0/64.0/gmPI/h9;
}

void SPHTension::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	attachAttribute(acceleration,std::string("ParticleAcceleration"));
	attachAttribute(density,std::string("ParticleDensity"));
	attachAttribute(plocality,std::string("ParticleLocalityGrid"));
}

int SPHTension::qlen()
{
	return 4;
}

void SPHTension::getq(double *q)
{
	q[0] = mass;
	q[1] = kernel_h;
	q[2] = tensionC;
	q[3] = normal_thresh;
}

void SPHTension::setq(double *q)
{
	mass = q[0];
	kernel_h = q[1];
	tensionC = q[2];
	normal_thresh = q[3];
}

void SPHTension::qname(char **qn)
{
	qn[0] = "mass";
	qn[1] = "kernel_h";
	qn[2] = "tension coefficient";
	qn[3] = "normal_threshhold";
}

gmVector3 SPHTension::dW_poly6(gmVector3 &r)
{
	double x=h2-dot(r,r);
	if (x>0)
		return -6.0*x*x*r*coeff;
	else
		return gmVector3(0,0,0);
}

double SPHTension::ddW_poly6(gmVector3 &r)
{
	double x=h2-dot(r, r);
	if (x>0)
		return -6.0*coeff*(3.0*x-4.0*dot(r,r) )*x;
	else
		return 0;
}

void SPHTension::applyForce()
{
	unsigned int i, k;
	int j;
	gmVector3 rij;  // vector between particle i and j
	gmVector3 normal;
	double kapa; // k

	std::vector<unsigned int> neighbors;
//	plocality->queryRadius=kernel_h;
	plocality->queryRadius=kernel_h*2;

  	for(i=0; i<ps->size(); i++)
	{
		normal = gmVector3(0.0, 0.0, 0.0);
		kapa = 0.0;

		neighbors.clear();
		plocality->getNeighbors(i, neighbors);
	
		for (k = 0; k < neighbors.size(); k++) 
		{
			j = neighbors[k];
//		for (j=0; j<ps->size(); j++)
//		{
			if (i==j) continue;   // not consider itself

			rij = position->getPosition(i) - position->getPosition(j);
			normal += (mass / density->den[j] ) * dW_poly6(rij); 
			kapa -= (mass / density->den[j] ) * ddW_poly6(rij);
		}

		if ( normal.length() > normal_thresh )
		{
			kapa /= normal.length();
			acceleration->acc[i]  += tensionC * kapa * normal / density->den[i];
		}
	}
}

void SPHTension::cleanup()
{
	if (plocality)
		plocality->update();
}
