/*
@file ParticleCreation.cpp
@author Wen Su
*/
#include "ParticleBoundingBox.h"

#include "ParticleCreation.h"
#include "ParticlePosition.h"


#ifdef _WIN32
#include <float.h>
#else
#include <math.h>
extern "C" int isnan(double);
#define _isnan isnan
#endif
		
REGISTER_PARTICLESTUFF(ParticleCreation,"Behavior:ParticleCreation");

ParticleCreation::ParticleCreation(Particles *ps, const std::string& name)
	: ParticleBehavior(ps, name)
{
	new PSParamInt(this,&limitSize,30,"minpop","min population",
		"Population threshold under which new particles are automatically generated.");

	new Attached<ParticleBoundingBox>(this,&bounds,"ParticleBoundingBox",
		"bounds","Bounding Box",
		"New particles will be created at uniform random locations in the region specified "
		"by this ParticleBounds attribute.");

	new Attached<ParticlePosition>(this,&position);
}

void ParticleCreation::cleanup()
{
	// add random pariticles;
	 gmVector3 &min=bounds->min;
	 gmVector3 &max=bounds->max;
	
	 if(_isnan(min[0]) || _isnan(min[1]) || _isnan(min[2]) ||
	   _isnan(max[0]) || _isnan(max[1]) || _isnan(max[2])){
		min = gmVector3(-10, -10, -10);
		max = gmVector3(10, 10, 10);
	}
	
	while(ps->size()<(unsigned int) limitSize)
	{
		ps->addParticle();
		position->setPosition(ps->size()-1,gmVector3(
			(max[0]-min[0])*(rand()/(double)RAND_MAX)+min[0],
			(max[1]-min[1])*(rand()/(double)RAND_MAX)+min[1],
			(max[2]-min[2])*(rand()/(double)RAND_MAX)+min[2]));
	}
}
