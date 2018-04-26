/**
 * Implementation of AdaptiveRepulsionData.
 * @file AdaptiveRepulsionData.cpp
 * @author Ed Bachta
 */

#include "libgm/gmConst.h"

#include "AdaptiveRepulsionData.h"

REGISTER_PARTICLESTUFF(AdaptiveRepulsionData,"Attribute:AdaptiveRepulsionData");

// Witkin-Heckbert recommended values.
#define WH_ALPHA 6.0

AdaptiveRepulsionData::AdaptiveRepulsionData(Particles* ps,const std::string& name)
: ParticleAttribute(ps, name)
{
	new PSParamDouble(this,&sigma_hat,0.25,
		"sigma_hat","desired repulsion radius",
		"Half of the desired distance between particles.");

	new PSParamDouble(this,&sigma_max,0.5,
		"sigma_max","max repulsion radius",
		"Half of the maximum distance between particles.");

	new PSParamDouble(this,&diameter,1.0,
		"diameter","surface diameter",
		"Estimate of the diameter of the surface.");

	new PSParamDouble(this,&alpha,6.0,
		"alpha","repusion amplitude",
		"Amplitude of the repulsion force. Proportional to particle motion response time, "
		"but if too large may cause stiffness");

	new PSParamDouble(this,&Ehat,0.8*alpha,
		"Ehat","desired energy",
		"Desired energy of the system. Too much energy leads to jittery particles. "
		"Too little causes the particle system to appear sluggish.");

	new PSParamDouble(this,&sdmul,3.0,
		"sdmul","repulsion clamp",
		"Particles beyond this factor of the radius are not repelled.");

	new PSParamDoublePerParticle(this,&r,"radius","particle radius",
		"width of Gaussian used to repel neighboring particles");

	new PSParamDoublePerParticle(this,&dr,"dradius","change in radius",
		"growth or reduction in particle radius");

	new PSParamDoublePerParticle(this,&D,"energy","particle energy",
		"energy accumulator for each particle");
}

/** Sets the particle system for which the
 * AdaptiveRepulsionData applies.
 * If new_ps is not NULL, then the lengths of internal arrays
 * of particle data are set to the size of the the particle
 * system.
 */
void AdaptiveRepulsionData::setParticleSystem(Particles *new_ps)
{
	ParticleAttribute::setParticleSystem(new_ps);
	if (ps)
	{
		r.resize(ps->size(), 1.0);
		dr.resize(ps->size(), 0.0);
		D.resize(ps->size(), 0.0);
	}
}



/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see ParticleAttribute::particleRemoved
 */
void AdaptiveRepulsionData::particleRemoved(unsigned int i)
{

	// Copy end element to ith position
	r[i]  = r.back();
	dr[i] = dr.back();
	D[i]  = D.back();

	// Pop off the end
	r.pop_back();
	dr.pop_back();
	D.pop_back();
}

/** 
 * Callback for particle addition.
 * @param i Index of the particle that has been added.
 */
void AdaptiveRepulsionData::particleAdded() 
{
	double rad = 1.0;

	r.push_back(rad);
	dr.push_back(0.0);
	D.push_back(0.0);
} 

void AdaptiveRepulsionData::clear()
{
	r.clear();
	dr.clear();
	D.clear();
}

void AdaptiveRepulsionData::integrate(double dt)
{
	for(unsigned int i=0;i<r.size();i++) 
	{
		r[i] += dr[i]*dt;
		// Seems to be more stable if I do not make all of them positive here
		r[i] = fabs(r[i]);
		//if (r[i] > sigma_max)
		//	r[i] = sigma_max;
	}
}
