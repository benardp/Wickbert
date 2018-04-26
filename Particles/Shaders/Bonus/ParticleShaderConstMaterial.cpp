#include "ConstantMaterial.h"
#include "ParticleMaterial.h"

REGISTER_PARTICLESTUFF(ConstantMaterial,"Shader:ConstantMaterial");

ConstantMaterial::ConstantMaterial(Particles *ps)
	:ParticleShader(ps,std::string("ConstantMaterial"))
{
	new PSParamgmVector4(this,&diffuseFront,gmVector4(1.0,0.0,0.0,1.0),
		"diff_front","Diffuse Front","Front facing diffuse color.");
	new PSParamgmVector4(this,&diffuseBack,gmVector4(0.0,0.0,1.0,1.0),
		"diff_back","Diffuse Back","Back facing diffuse color.");

	new Attached<ParticleMaterial>(this,&material);
}

void ParticleShaderConstMaterial::drawShape(int i)
{
	if (material) {
		material->setDiffuseFront(i,diffuseFront);
		material->setDiffuseBack(i,diffuseBack);
		material->color[i]=color;
	}
}
