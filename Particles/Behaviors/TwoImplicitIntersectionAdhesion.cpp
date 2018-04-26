#include "TwoImplicitIntersectionAdhesion.h"
#include "ParticleAge.h"


TwoImplicitIntersectionAdhesion::TwoImplicitIntersectionAdhesion(Particles *ps, const std::string & name)
: ParticleBehavior(ps, name)
{
	new Attached<ParticleAge>(this,&pAge);
}

void TwoImplicitIntersectionAdhesion::attachAttributes()
{
	ParticleBehavior::attachAttributes();
}

