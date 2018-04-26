#include "UseMaterial.h"
#include "ParticleMaterial.h"

REGISTER_PARTICLESTUFF(UseMaterial,"Shader:UseMaterial");

UseMaterial::UseMaterial(Particles *ps)
	:ParticleShader(ps,std::string("UseMaterial"))
{
	new Attached<ParticleMaterial>(this,&material);
}

void UseMaterial::drawParticle(int i)
{
	material->setMaterial(i);
}


