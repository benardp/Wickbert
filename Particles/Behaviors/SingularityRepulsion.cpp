 /*
@file SingularityRepulsion.cpp
@author Wen Su
*/

#include "SingularityRepulsion.h"

#include "pstools.h"
#include <math.h>
#include "Particles.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
//#include "ImplicitInterrogator.h"
#include "ParticleRepulsion.h"
#include "AdaptiveRepulsionData.h"

REGISTER_PARTICLESTUFF(SingularityRepulsion,"Behavior:SingularityRepulsion");

SingularityRepulsion::SingularityRepulsion(Particles *ps /*, Particles *t=NULL , ParticleRepulsion *r=NULL*/)
	: ParticleBehavior(ps, std::string("SingularityRepulsion"))
{
	//target=NULL;
	// find repulsion in target system
	repulsion=NULL;
	targetRepelConst=10;
	selfRepelConst=0.0005;	// selfRepelsion

	new PSParamString(this,&target,"floaters","target","target particles",
					  "Name of Particles ...");
	
	new Attached<ParticlePosition>(this,&pos);
	new Attached<ParticleVelocity>(this,&vel);
}

void SingularityRepulsion::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	if (!ps) return;
	target_p = ps->particleSystem->findParticles(target);
	//impInt=target_p->getBehavior<ImplicitInterrogator>(std::string("ImplicitInterrogator"));
}


int SingularityRepulsion::qlen()
{
	return 2;
}

void SingularityRepulsion::getq(double *q)
{
	q[0] = targetRepelConst;
	q[1] = selfRepelConst;
}

void SingularityRepulsion::setq(double *q)
{
	targetRepelConst = q[0];
	selfRepelConst = q[1];
}

void SingularityRepulsion::qname(char **qn)
{
	qn[0] = "Target Repulsion constant";
	qn[1] = "Self Repulsion constant";
}

void SingularityRepulsion::applyForce()
{
	if (target_p==NULL || repulsion==NULL)
		return;

	repelTarget();
	repelSelf();
}

/// this should be exactly the same as particle repulsion
void SingularityRepulsion::repelTarget()
{
	double rij2; // distance squared between particle and one of its neighbors
	double sig2; // particle's repulsion radius squared
	double Eij;	// Energy term between particle and one of its neighbors
	gmVector3 rij; // Vector between particles i and j.
	
	ParticlePosition *targetPosition=target_p->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
	ParticleVelocity *targetVelocity=target_p->getAttribute<ParticleVelocity>(std::string("ParticleVelocity"));
	AdaptiveRepulsionData* rep_data=target_p->getAttribute<AdaptiveRepulsionData>(std::string("AdaptiveRepulsionData"));
	// for all singular points in this system
	for (unsigned int i = 0; i < ps->size(); i++)
	{
		// Apply the force based on all neighbors
		for (unsigned int j = 0; j < target_p->size(); j++)
		{
			// self - target
			rij = pos->getPosition(i) - targetPosition->getPosition(j);
			rij2 = rij.lengthSquared();
			sig2 = rep_data->r[i];
			sig2 = sig2 * sig2;
			
			// Compute energy contribution
			Eij = rep_data->alpha*fastExp(-0.5*rij2/sig2);

			// Scale Eij by distance so it goes to zero when |rij|
			Eij *= 1.0 - rij2/(sig2*rep_data->sdmul*rep_data->sdmul);
			
			// Add energy into particle j's velocity (second half of (9))
			targetVelocity->v[j] -= targetRepelConst*(rep_data->r[j]*rep_data->r[j]/sig2)*rij*Eij;
		}
	}
}

/// this is a simpler repulsion because singular points has less population
void SingularityRepulsion::repelSelf()
{
	double rij2; // distance squared between particle and one of its neighbors
	double sig2=1; // particle's repulsion radius squared
	gmVector3 rij; // Vector between particles i and j.
	gmVector3 force;

	// for all singular points in this system
	for (unsigned int i = 0; i < ps->size(); i++)
	{
		// no repeatition
		for (unsigned int j = 0; j < i; j++)
		{
			rij=pos->getPosition(i)-pos->getPosition(j);
			rij2=rij.lengthSquared();
			if (rij2<sig2)
			{
				force=selfRepelConst*rij/rij2;
				vel->v[i]+=force;
				vel->v[j]-=force;
			}
		}
	}
}

