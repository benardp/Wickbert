/**
 * @file TopologyInterrogator.cpp
 * @date 2003-09-07
 * @author Wen Su
 */

#include "TopologyInterrogator.h"
#include "ImplicitInterrogator.h"
#include "ParticleMesh.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"

REGISTER_PARTICLESTUFF(TopologyInterrogator,"Behavior:TopologyInterrogator");

TopologyInterrogator::TopologyInterrogator(Particles *ps)
	: ParticleBehavior(ps, std::string("TopologyInterrogator"))
{
	im=NULL;
	reStart=true;
	startUpSize=1.0;
}

void TopologyInterrogator::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	attachAttribute(im,std::string("ImplicitInterrogator"));
}

void TopologyInterrogator::applyConstraint()
{
	gmVector3 xi,gradient;
	double scale;
	// push particles to surface;
	for(unsigned int i=0;i<ps->size();++i)
	{
		xi = position->getPosition(i);
		gradient=im->getImplicit()->grad(xi);
		if (gradient.lengthSquared() < 0.001)
			continue;
		scale = 0;//(dot(gradient, ps->v[i])+15*im->getImplicit()->proc(xi))/gradient.lengthSquared();
		velocity->v[i] -= scale*gradient;
	}
}

void TopologyInterrogator::applyForce()
{
}

void TopologyInterrogator::cleanup()
{
	ParticleMesh *pm=dynamic_cast<ParticleMesh *>(position);

	if (pm==NULL)
		return;
	/// todo: this really should be called from Topology change when ready
	if (reStart)
	{
		reStart=false;
		ps->removeAll();
		pm->createDiamond(startUpSize);
	}
}

int TopologyInterrogator::qlen()
{
	return 2;
}
void TopologyInterrogator::getq(double *q)
{
	q[0] = startUpSize;
	q[1] = reStart;
}
void TopologyInterrogator::setq(double *q)
{
	startUpSize = q[0];
	reStart = q[1];
}
void TopologyInterrogator::qname(char **qn)
{
	qn[0] = "start up size";
	qn[1] = "restart";
}
