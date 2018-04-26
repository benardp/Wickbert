/**
* @file SPHViscosity.cpp
* @author Yuan Zhou
*/

#include "SPHViscosity.h"

#define mass_default 1.0 
#define kernel_h_default 1.0 //0.11 // 1.0
#define viscosity_default 1.0 // 2.0 // 0.1  // 0.5 // 100   

REGISTER_PARTICLESTUFF(SPHViscosity,"Behavior:SPHViscosity");

SPHViscosity::SPHViscosity(Particles *ps) 
: ParticleBehavior(ps, std::string("SPHViscosity") )
{
    density = NULL;
	acceleration = NULL;
	mass = mass_default;
	kernel_h = kernel_h_default;
	viscosity = viscosity_default;
	plocality=NULL;

	h2 = kernel_h*kernel_h;
	h6 = h2*h2*h2;
	coeff = 45.0/gmPI/h6;
}

void SPHViscosity::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	attachAttribute(acceleration,std::string("ParticleAcceleration"));
	attachAttribute(density,std::string("ParticleDensity"));
	attachAttribute(plocality,std::string("ParticleLocalityGrid"));
}

int SPHViscosity::qlen()
{
	return 3;
}

void SPHViscosity::getq(double *q)
{
	q[0] = mass;
	q[1] = kernel_h;
	q[2] = viscosity;
}

void SPHViscosity::setq(double *q)
{
	mass = q[0];
	kernel_h = q[1];
	viscosity = q[2];
}

void SPHViscosity::qname(char **qn)
{
	qn[0] = "mass";
	qn[1] = "kernel_h";
	qn[2] = "viscosity";
}

double SPHViscosity::ddW_visco(gmVector3 &r)
{
	double r2 = dot(r,r);
	if (r2 < h2)
		return coeff*( kernel_h - sqrt(r2) );
	else
		return 0;
}


void SPHViscosity::applyForce()
{
	unsigned int i, k;
	int j;
	gmVector3 rij;  // vector between particle i and j
	gmVector3  total;

	std::vector<unsigned int> neighbors;
//	plocality->queryRadius=kernel_h;
	plocality->queryRadius=kernel_h*2;

    for(i=0; i<ps->size(); i++)
	{
		total  = gmVector3(0.0, 0.0, 0.0);

		neighbors.clear();
		plocality->getNeighbors(i, neighbors);

		for (k = 0; k < neighbors.size(); k++) 
		{
			j = neighbors[k];
//		for (j=0; j<ps->size(); j++)
//		{
			if ( j==i) continue;
			rij = position->getPosition(i) - position->getPosition(j);
			total += ( mass * (velocity->v[j] - velocity->v[i]) / density->den[j] ) * ddW_visco(rij);
		}
		acceleration->acc[i] += viscosity*total / density->den[i];  
	}
}

void SPHViscosity::cleanup()
{
	if (plocality)
		plocality->update();
}
