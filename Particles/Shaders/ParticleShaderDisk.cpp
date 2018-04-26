#include "ParticleShaderDisk.h"


REGISTER_PARTICLESTUFF(ParticleShaderDisk,"Shader:ParticleShaderDisk");

ParticleShaderDisk::ParticleShaderDisk(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderDisk"))
{
	quad=gluNewQuadric();

	new PSParamString(this,&radius_source,"AdaptiveRepulsionData:radius","radattr","radius source",
		"Source of the radius of the form Attribute:param");
	new PSParamDouble(this,&scale,0.5,"scale","scale factor",
		"Amount to scale the radius");
	new PSParamInt(this,&sides,10,"sides","polygon sides",
		"# of sides of the disk approximating polygons");

	radius_attr = NULL;
	radius_data = NULL;
}

void ParticleShaderDisk::attachAttributes()
{
	PSParamPerParticle *radius_param = NULL;

	ParticleShader::attachAttributes();

	std::string attr_name = radius_source.substr(0,radius_source.find(':'));
	std::string param_name = radius_source.substr(radius_source.find(':')+1);
	if (radius_attr = ps->getAttribute<ParticleAttribute>(attr_name))
		radius_param = radius_attr->perparticle.findparam(param_name);
	if (radius_param)
		radius_data = (DoubleVector *)radius_param->ref();
	else
		radius_attr = NULL;

#if 0
	if (radius_attr) {
		int i;
		for (i = 0; i < radius_attr->perparticle.size(); i++) {
			if (radius_attr->perparticle[i]->shortname == param_name) {
				radius_data = (DoubleVector *)radius_attr->perparticle[i]->ref();
				break;
			}
		}
		// If param not found, toss attribute too
		if (i == radius_attr->perparticle.size())
			radius_attr = NULL;
	}
#endif
}

void ParticleShaderDisk::drawParticle(int i)
{
	if (radius_data)
		radius = (*radius_data)[i];
	else
		radius = 1.0;

	gluDisk(quad,0,radius*scale,sides,sides);
}

ParticleShaderDisk::~ParticleShaderDisk()
{
	gluDeleteQuadric(quad);
}
