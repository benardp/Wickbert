//@file SPHConstraint.cpp
//@author Yuan Zhou

#include "SPHConstraint.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"

#define safe 0.0 // 1e-8

REGISTER_PARTICLESTUFF(SPHConstraint,"Behavior:SPHConstraint");

SPHConstraint::SPHConstraint(Particles* ps, const std::string& name)
:ParticleBehavior(ps, name)
{
	bbox=NULL;
}

char *SPHConstraint::tip()
{
	return "Restricts particle positions to a bounding box.";
}

void SPHConstraint::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	attachAttribute(bbox,std::string("ParticleBoundingBox"));
}

void SPHConstraint::applyConstraint()
{
	if (bbox==NULL)
		return;
	gmVector3 &min=bbox->min;
	gmVector3 &max=bbox->max;

	for (unsigned int i=0;i<ps->size();i++)
	{		
		if (position->getPosition(i)[2]<=min[2])  // bottom
		{			
			velocity->v[i][2] = -velocity->v[i][2];
			position->x[i][2] = min[2]+safe;
		}
		if (position->getPosition(i)[2] > bbox->max[2]) // top
		{
			velocity->v[i][2] = -velocity->v[i][2];
			position->x[i][2] = max[2]-safe;
		}
		
		if (position->getPosition(i)[1] <=min[1])  // left
		{ 
			velocity->v[i][1] = - velocity->v[i][1];
			position->x[i][1] = min[1] +safe;
		}
		if (position->getPosition(i)[1] >= max[1]) // right
		{
			velocity->v[i][1] = - velocity->v[i][1];
			position->x[i][1] = max[1]-safe;
		}

		if (position->getPosition(i)[0] <= min[0])   // rear
		{
			velocity->v[i][0] = - velocity->v[i][0];
			position->x[i][0] = min[0]+safe;  
		}
		if (position->getPosition(i)[0] >= max[0])  // front
		{
			velocity->v[i][0] = - velocity->v[i][0];
			position->x[i][0] = max[0] - safe;
		}
	} 
}
