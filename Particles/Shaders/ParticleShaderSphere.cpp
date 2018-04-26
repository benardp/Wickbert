#include "ParticleShaderSphere.h"
#include "AdaptiveRepulsionData.h"

REGISTER_PARTICLESTUFF(ParticleShaderSphere,"Shader:ParticleShaderSphere");

ParticleShaderSphere::ParticleShaderSphere(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderSphere"))
{
	quad=gluNewQuadric();

	new PSParamDouble(this,&radius,0.1,"r","Radius","Radius of sphere.");
	new PSParamInt(this,&sides,10,"n","Sides","# of sides with which to polygonize sphere.");
}

ParticleShaderSphere::~ParticleShaderSphere()
{
	gluDeleteQuadric(quad);
}


void ParticleShaderSphere::drawParticle(int i)
{
	gluSphere(quad,radius,sides,sides);
}
