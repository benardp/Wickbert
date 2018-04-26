#include "UseMaterialSwitch.h"
#include "ParticleMaterial.h"
#include "ParticleScalar.h"

REGISTER_PARTICLESTUFF(UseMaterialSwitch,"Shader:UseMaterialSwitch");

UseMaterialSwitch::UseMaterialSwitch(Particles *ps)
	:ParticleShader(ps,std::string("UseMaterialSwitch"))
{
	//	new PSParamString(this,&radius,0.1,"r","Radius","Radius of sphere.");
	new Attached<ParticleScalar>(this,&switchScalar);
	new Attached<ParticleMaterial>(this,&materialNegative);
	new Attached<ParticleMaterial>(this,&materialZero);
	new Attached<ParticleMaterial>(this,&materialPositive);
}

void UseMaterialSwitch::drawParticle(int i)
{
	double val = switchScalar->getScalar(i);
	std::string a = materialZero->name;
	if(val > 0)
		materialPositive->setMaterial(i);
	else if(val < 0)
		materialNegative->setMaterial(i);
	else
		materialZero->setMaterial(i);
}
