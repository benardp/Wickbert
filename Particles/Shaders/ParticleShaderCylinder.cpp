#include "ParticleShaderCylinder.h"

REGISTER_PARTICLESTUFF(ParticleShaderCylinder,"Shader:ParticleShaderCylinder");

ParticleShaderCylinder::ParticleShaderCylinder(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderCylinder"))
{
	quad=gluNewQuadric();
	base=0.1f;
	height=0.3f;
	top=0;
	sides=10;
	capped = 1;
}

ParticleShaderCylinder::~ParticleShaderCylinder()
{
	gluDeleteQuadric(quad);
}

void ParticleShaderCylinder::drawParticle(int i)
{
	float radius=1.0f;
	
	gluCylinder(quad,radius*base,radius*top,height,sides,1);
	if (capped) {
		gluDisk(quad,0.0f,radius*base,sides,1);
		glTranslatef(0.0,0.0,height);
		gluDisk(quad,0.0,radius*top,sides,1);
	}
}
