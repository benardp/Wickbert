#include "SingularityAdhesion.h"
#include "ImplicitInterrogator.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "Surface/Implicit/Implicit.h"


REGISTER_PARTICLESTUFF(SingularityAdhesion,"Behavior:SingularityAdhesion");

/// default constructor
SingularityAdhesion::SingularityAdhesion(Particles* ps, const std::string& name)
	: ParticleBehavior(ps, name)
{
	new PSParamDouble(this,&strength,10.0,"strength","force strength",
		"Scalar multplied times gradient of gradient magnitude.");
	new PSParamDouble(this,&clamp,0.2,"clamp","force clamp",
		"Maximum magnitude of gradient of gradient magnitude.");

	new Attached<ImplicitInterrogator>(this,&impInt);
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
}

// this approaches a singularity point
void SingularityAdhesion::applyForce()
{
	if (!impInt) return;

	// move towards the opposite of the gradient direction.
	for (unsigned int i = 0; i < ps->size(); i++)
	{
		gmVector3 gradGradN2=2.0*impInt->hess(i)*impInt->grad(i);
		gmVector3 dir = gradGradN2;//impInt->getImplicit()->gradGradN2(position->getPosition(i));
		dir *= -strength;
		double length = dir.length();
		if (clamp < length)
			dir *= clamp/length;
		velocity->v[i] += dir;
	}
}
