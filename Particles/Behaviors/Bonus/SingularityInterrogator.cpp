#include "SingularityInterrogator.h"
#include "ImplicitInterrogator.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "Surface/Implicit/Implicit.h"


REGISTER_PARTICLESTUFF(SingularityInterrogator,"Behavior:SingularityInterrogator");

/// default constructor
SingularityInterrogator::SingularityInterrogator(Particles* ps, const std::string& name)
	: ParticleBehavior(ps, name)
{
	gradientConst=10;
	feedbackConst=15;
	initIntervalSize=1;
	maxLength=0.2;
	impInt=NULL;
}

void SingularityInterrogator::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	attachAttribute(impInt,std::string("ImplicitInterrogator"));
}

int SingularityInterrogator::qlen()
{
	return 3;
}

void SingularityInterrogator::getq(double *q)
{
	q[0] = gradientConst;
	q[1] = feedbackConst;
	q[2] = maxLength;
}

void SingularityInterrogator::setq(double *q)
{
	gradientConst = q[0];
	feedbackConst = q[1];
	maxLength = q[2];
}

void SingularityInterrogator::qname(char **qn)
{
	qn[0] = "gradient scaling constant";
	qn[1] = "surface feedback scaling constant";
	qn[2] = "maximum length for movement";
}

// this approaches a singularity point
void SingularityInterrogator::applyForce()
{
	if (!((impInt) && (impInt->getImplicit())))
		return;
	// move towards the opposite of the gradient direction.
	for (unsigned int i = 0; i < ps->size(); i++)
	{
		gmVector3 dir=impInt->getImplicit()->gradGradN2(position->getPosition(i));
		dir=-gradientConst*dir;
		double length=dir.length();
		if (maxLength<length)
		{
			dir.normalize();
			dir=maxLength*dir;
		}
		velocity->v[i]+=dir;
	}
}
