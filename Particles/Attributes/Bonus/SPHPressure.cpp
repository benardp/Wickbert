/**
* @file SPHPressure.cpp
* @author Yuan Zhou
*/

#include "SPHPressure.h"

#define mass_default 1.0 //1.0 
#define gasC_default   50.0 // 25.0 // 1000.0    // gas constant for pressure
#define pressure_kernel_default  1.0 // 0.11 
#define density_rest_default 0.0   // 1000.0 

REGISTER_PARTICLESTUFF(SPHPressure,"Behavior:SPHPressure");

SPHPressure::SPHPressure(Particles *ps) 
: ParticleBehavior(ps, std::string("SPHPressure") )
{
	density = NULL;
	acceleration = NULL;
	mass = mass_default;
	gasC = gasC_default;
	pressure_kernel = pressure_kernel_default;
	plocality=NULL;
	density_rest = density_rest_default;

	h2 = pressure_kernel * pressure_kernel;
	h6 = h2 * h2 * h2;
	coeff = 45.0/gmPI/h6;
}

void SPHPressure::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	attachAttribute(acceleration,std::string("ParticleAcceleration"));
	attachAttribute(density,std::string("ParticleDensity"));
	attachAttribute(plocality,std::string("ParticleLocalityGrid"));
}

int SPHPressure::qlen()
{
	return 4;
}

void SPHPressure::getq(double *q)
{
	q[0] = mass;
	q[1] = gasC;
	q[2] = pressure_kernel;
	q[3] = density_rest;

}

void SPHPressure::setq(double *q)
{
	mass = q[0];
	gasC = q[1];
	pressure_kernel = q[2];
	density_rest = q[3];
}

void SPHPressure::qname(char **qn)
{
	qn[0] = "mass";
	qn[1] = "gasC";
	qn[2] = "pressure_kernel_h";
	qn[3] = "density_rest";
}

gmVector3 SPHPressure::dW_spiky(gmVector3 &r)
{
	double r2 = dot(r, r);
	if ( sqrt(r2) < 1e-6*pressure_kernel || r2 > h2 )
		return gmVector3(0,0,0);
	double x = sqrt(r2) / pressure_kernel;
	return (-coeff*pressure_kernel*(1-x)*(1-x)/x)*r;
}

void SPHPressure::applyForce()
{
	unsigned int i, k;
	int j;

	gmVector3 rij;  // vector between particle i and j
	double sqr_ph = pressure_kernel * pressure_kernel;  // h2
	double total;
	gmVector3  pressure;

	std::vector<unsigned int> neighbors;

//	plocality->queryRadius=pressure_kernel;
	plocality->queryRadius=pressure_kernel*2;

 	// compute pressure gradient force for all particles
	for(i=0; i<ps->size(); i++)
	{
		pressure = gmVector3(0.0, 0.0, 0.0);
		double pre;
		neighbors.clear();
		plocality->getNeighbors(i, neighbors);
		std::cout << "#neighbor: " << neighbors.size() << std::endl;

		for (k = 0; k < neighbors.size(); k++) 
		{
			j = neighbors[k];

//		for(j=0; j<ps->size(); j++) 
//		{
			if (i==j) continue;  // shouldn't consider itself
			rij = position->getPosition(i) - position->getPosition(j);
			pre = mass * gasC * (density->den[i] + density->den[j] - 2.0*density_rest) / (2.0*density->den[j]);

			pressure -= pre * dW_spiky(rij);
		}
		acceleration->acc[i] += pressure / density->den[i];  
	}
}

void SPHPressure::cleanup()
{
	if (plocality)
		plocality->update();
}