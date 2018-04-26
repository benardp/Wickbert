/**
 * Implementation of ParticleFate.
 * @file ParticleFate.cpp 
 * @date November 28, 2001 
 * @author Ed Bachta
 */

#include "ParticleFate.h"
#include "pstools.h"
#include "AdaptiveRepulsionData.h"
#include "ParticleNormal.h"
#include "ParticleBoundingBox.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "ImplicitInterrogator.h"
#include "ParticleAge.h"

namespace{
	double max(double a,double b)
	{
		return (a>b)?a:b;
	}	
}

#define INITIAL_DIAMETER 1.0

REGISTER_PARTICLESTUFF(ParticleFate,"Behavior:ParticleFate");

/**
 * Creates a ParticleFate object. First it looks for the required 
 * attributes in the particle system, and then it sets default 
 * constants.
 */
ParticleFate::ParticleFate(Particles *ps, const std::string &name) 
	: ParticleBehavior(ps, name)
{
	new PSParamDouble(this,&gamma,4.0,"gamma","Equilibrium rate",
		"Equilibrium velocity as a factor of repulsion radius");
	new PSParamDouble(this,&nu,0.2,"nu","fission threshold",
		"Particle divides when its energy exceeds this fraction of desired energy");
	new PSParamDouble(this,&delta,0.7,"delta","death threshold",
		"Particle dies when its radius falls below this fraction of the desired repulsion radius");
	new PSParamInt(this,&population,100,"cap","max population","Upper limit on the population");
	new PSParamDouble(this,&fateval,0.1,
		"fateval","Fate Threshold","Fission/death only occur when absolute particle value within this threshold.");

	new Attached<AdaptiveRepulsionData>(this,&rep_data);
	new Attached<ParticleNormal>(this,&p_orient);
	new Attached<ImplicitInterrogator>(this,&impint,"ImplicitInterrogator",
		"impint","implicit interrogator","Fission/death only occur near implicit surface.");
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
	new Attached<ParticleAge>(this,&pAge);
}

void ParticleFate::attachAttributes()
{
	ParticleBehavior::attachAttributes();
}

/** 
* Allows for an awareness of surface rep_data->diameter.
* @param d New surface rep_data->diameter.
* @returns New desired readius.
*/
double ParticleFate::setSurfaceDiameter(double d)
{
	// rep_data->sigma_hat and rep_data->sigma_max are updated as suggested by Witkin and Heckbert
	if (d > 0.0) {
		rep_data->diameter = d;
		rep_data->sigma_hat = gmMin(rep_data->sigma_hat, rep_data->diameter/4.0);
		rep_data->sigma_max = gmMax(rep_data->diameter/2.0, 1.5*rep_data->sigma_hat);
	} 
	return rep_data->sigma_hat;
}

/**
* Allow the user to specify the desired radius. This will modify
* the number and size of particles on the surface.
* @param r Desired radius.
* @note The resulting radius is bounded by the surface rep_data->diameter.
*/
void ParticleFate::setDesiredRadius(double &r)
{
	// rep_data->sigma_hat and rep_data->sigma_max are updated as suggested by Witkin and Heckbert
	if (r<=0.0)
	{
		r = rep_data->sigma_hat;
		return;
	}
	else
	{
		r = gmMin(r, rep_data->diameter/4.0);
		rep_data->sigma_hat = r;
		rep_data->sigma_max = gmMax(rep_data->diameter/2.0, 1.5*rep_data->sigma_hat);
	}
}

/**
* We determine particle fate during the cleanup step.
*/
void ParticleFate::cleanup()
{
	Implicit *imp;

	if (impint)
		imp = impint->getImplicit();
	else
		imp = NULL;

	// Iterate backward to preserve indexing while removing particles
	for (int i=ps->size()-1; i>=0; i--) 
	{
		if (imp) {
			double f = impint->proc(i);//imp->proc(position->x[i]);
			if (fabs(f) > fateval)
				continue;
		}

		// check if radius is negative
		if (rep_data->r[i] < gmEPSILON)
			ps->removeParticle(i); //realize that this also updates the implicit interrogator

		/**
		 * Determine whether or not particle is to die:
		 * 1) Particle velocity is small compared to nominal repulsion radius
		 * 2) Repulsion radius is smaller than minimum repulsion radius
		 * 3) Biased random test succeds
		 * OR
		 * 4) Particle falls outside of the bounding box
		 */
		if ((velocity->v[i].length() < gamma*rep_data->r[i]) &&
			(rep_data->r[i] < delta*rep_data->sigma_hat) &&
			(psRandom() > rep_data->r[i] / (delta*rep_data->sigma_hat)))
			ps->removeParticle(i);

		if(pAge->t[i] > pAge->maxLifeTime) ps->removeParticle(i);
		
		// this is not your job
		//if (p_bounds)
		//	p_bounds->keepInBounds(position->getPosition(i));
	}
	
	// Iterate forward to new size and split particles
	// ps size may change, only allow one split per particle
	unsigned int oldSize=ps->size();
	for (unsigned int i=0; i<oldSize; i++) 
	{
		// max population
		if (oldSize>(unsigned int)population || i>(unsigned int)population)
			break;

		if (imp) {
			double f = impint->proc(i);//imp->proc(position->x[i]);
			if (fabs(f) > fateval)
				continue;
		}

		// Make sure that particles with no direction to grow
		// don't spawn new particles, ie. no active particles
		// for base cases like "Difference". 
		if ((p_orient) && (p_orient->getGradMag2(i) < gmEPSILON))
			continue;

		/*
		 * Determine whether or not to split:
		 * (Condition 1:  particle is near equilibrium) AND
		 * (Condition 2a: repulsion radius is huge OR
		 *  Condition 2b: particle is adequately energized and repulsion radius
		 *				 is higher than nominal)
		 */
		if (
			(velocity->v[i].length() < gamma*rep_data->r[i]) 
			&&
			(
				(rep_data->r[i] > rep_data->sigma_max) || 
				(
					(rep_data->D[i] > nu*rep_data->Ehat) && 
					(rep_data->r[i] > rep_data->sigma_hat)
				)
			) 
		) 
		{
			
			// Caclulate perturbation
			double pert = 0.5*rep_data->r[i];
			
			// Generate tangent vectors orthogonal to the normal
			gmVector3 tgt1 = randOrtho(p_orient->getNormal(i));
			tgt1.normalize();
			
			gmVector3 tgt2 = cross(p_orient->getNormal(i), tgt1);
			tgt2.normalize();
			
			gmVector3 disp;
			disp = tgt1 * (2.0 * pert * psRandom() - pert) +
				tgt2 * (2.0 * pert * psRandom() - pert);
			
			gmVector3 newLoc = position->getPosition(i) + disp;
			// this is not your job
			//if (p_bounds)
			//	p_bounds->keepInBounds(newLoc);
			
			// The new particle is added
			int j=ps->size();
			ps->addParticle();
			position->setPosition(j,newLoc);
			
			// Reduce the repulsion radii for both particles
			rep_data->r[i] /= sqrt(2.0);
			// The following line should not be necessary, but
			// you can uncomment if particles become ill-behaved.
			// rep_data->r[i] = fabs(rep_data->r[i]);
			rep_data->r[j] = rep_data->r[i];	
		}	
	}	
}
