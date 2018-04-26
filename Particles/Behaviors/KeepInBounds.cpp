//@file KeepInBounds.cpp
//@date 2004-02-19
//@author Wen Su

#include "KeepInBounds.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "ParticleBoundingBox.h"

REGISTER_PARTICLESTUFF(KeepInBounds,"Behavior:KeepInBounds");

KeepInBounds::KeepInBounds(Particles* ps, const std::string& name)
:ParticleBehavior(ps, name)
{
	new Attached<ParticleBoundingBox>(this,&bbox);
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
}

void KeepInBounds::applyConstraint()
{
	if (bbox==NULL)
		return;
	gmVector3 &min=bbox->min;
	gmVector3 &max=bbox->max;

	for (unsigned int i=0;i<ps->size();i++)
	{		
		if (position->getPosition(i)[2] > bbox->max[2])
		{			
			// Top constraint
			velocity->v[i][0] = 0.0;
			velocity->v[i][1] = 0.0;
			velocity->v[i][2] = (velocity->v[i][2] > 0.0) ? 0.0 : velocity->v[i][2];			
		}
		else if (position->getPosition(i)[2] < min[2])
		{			
			// Bottom constraint
			velocity->v[i][0] = 0.0;
			velocity->v[i][1] = 0.0;
			velocity->v[i][2] = (velocity->v[i][2] < 0.0) ? 0.0 : velocity->v[i][2];			
		}
		else if (position->getPosition(i)[1] > max[1])
		{
			// Right constraint
			velocity->v[i][0] = 0.0;
			velocity->v[i][2] = 0.0;
			velocity->v[i][1] = (velocity->v[i][1] > 0.0) ? 0.0 : velocity->v[i][1];			
		}
		else if (position->getPosition(i)[1] < min[1])
		{			
			// Left constraint
			velocity->v[i][0] = 0.0;
			velocity->v[i][2] = 0.0;
			velocity->v[i][1] = (velocity->v[i][1] < 0.0) ? 0.0 : velocity->v[i][1];			
		}
		else if (position->getPosition(i)[0] > max[0])
		{			
			// Front constraint
			velocity->v[i][1] = 0.0;
			velocity->v[i][2] = 0.0;
			velocity->v[i][0] = (velocity->v[i][0] > 0.0) ? 0.0 : velocity->v[i][0];			
		}
		else if (position->getPosition(i)[0] < min[0])
		{			
			// Rear constraint
			velocity->v[i][1] = 0.0;
			velocity->v[i][2] = 0.0;
			velocity->v[i][0] = (velocity->v[i][0] < 0.0) ? 0.0 : velocity->v[i][0];
			
		} // end bounds checking
	} // end particle loop
}
