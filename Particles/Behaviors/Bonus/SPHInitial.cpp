/**
* @file SPHInitial.cpp
* @author Yuan Zhou
*/

#include "SPHInitial.h"

#define mass_default 1.0 
#define density_kernel_default 1.0 // 0.11 

REGISTER_PARTICLESTUFF(SPHInitial,"Behavior:SPHInitial");

SPHInitial::SPHInitial(Particles *ps) 
: ParticleBehavior(ps, std::string("SPHInitial") )
{
	density = NULL;
	acceleration = NULL;
	mass = mass_default;
	density_kernel = density_kernel_default;
	plocality=NULL;

	h2 = density_kernel * density_kernel;  // h^2
	h9 = h2*h2*h2*h2*density_kernel; // h^9
	coeff = 315.0/64.0/gmPI/h9;
}

void SPHInitial::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	attachAttribute(acceleration,std::string("ParticleAcceleration"));
	attachAttribute(density,std::string("ParticleDensity"));
	attachAttribute(plocality,std::string("ParticleLocalityGrid"));
}

int SPHInitial::qlen()
{
	return 2;
}

void SPHInitial::getq(double *q)
{
	q[0] = mass;
	q[1] = density_kernel;
}

void SPHInitial::setq(double *q)
{
	mass = q[0];
	density_kernel = q[1];
}

void SPHInitial::qname(char **qn)
{
	qn[0] = "mass";
	qn[1] = "density_kernel_h";
}

double SPHInitial::W_poly6(const gmVector3 &r)
{
	double x=h2-dot(r, r);
	if (x>0)
		return coeff*x*x*x;
	else
		return 0;
}

void SPHInitial::applyForce()
{
	unsigned int i, k;
	int j;

	gmVector3 rij;  // vector between particle i and j
	double total;
	
	std::vector<unsigned int> neighbors;
//	plocality->queryRadius=density_kernel;
	plocality->queryRadius=density_kernel*2;

    for(i=0; i<ps->size(); i++)
	{
		total = W_poly6( gmVector3(0,0,0) );  // density from itself
		neighbors.clear();
		plocality->getNeighbors(i, neighbors);
	//	printf("# neighbors: %d\n", neighbors.size());

		for (k = 0; k < neighbors.size(); k++) 
		{
			j = neighbors[k];
//		for (j=0; j<ps->size(); j++)
//		{
			if (i == j) continue;
			rij = position->getPosition(i) - position->getPosition(j);
			total += W_poly6(rij);
		}
		density->den[i] = mass*total;  // assume all particles have the same mass
		acceleration->acc[i] = gmVector3(0,0,0);
	}
}

void SPHInitial::cleanup()
{
	if (plocality)
		plocality->update();
}