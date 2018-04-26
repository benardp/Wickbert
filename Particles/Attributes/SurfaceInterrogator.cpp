#include "SurfaceInterrogator.h"
#include "Surface/Surface.h"
#include "ParticleSystem.h"

REGISTER_PARTICLESTUFF(SurfaceInterrogator,"Attribute:SurfaceInterrogator");

SurfaceInterrogator::SurfaceInterrogator(Particles* ps, Surface* s, const std::string& name)
	: ParticleAttribute(ps, name)
{
	//index=0;
	surface=s;

	new PSParamString(this,&surfname,"<empty>",
		"surface","target surface","The surface to interrogate.");
}

void SurfaceInterrogator::attachAttributes()
{
	surface = ps->particleSystem->surfaces->withName(surfname);
}

