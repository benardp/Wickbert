#include "ParticleShaderOriented.h"
#include "AdaptiveRepulsionData.h"

REGISTER_PARTICLESTUFF(ParticleShaderOriented,"Shader:ParticleShaderOriented");

ParticleShaderOriented::ParticleShaderOriented(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderOriented"))
{
	quad = gluNewQuadric();
	radius=0.1;
	scale=0.5;
	sides=10;
	imp_int = 0;
}

void ParticleShaderOriented::attachAttributes()
{
	ParticleShader::attachAttributes();
	attachAttribute(ARData,std::string("AdaptiveRepulsionData"));
	attachAttribute(imp_int,std::string("ImplicitInterrogator"));
}

int ParticleShaderOriented::qlen()
{
	return 2;
}

void ParticleShaderOriented::getq(double *q)
{
	q[0] = sides;
	q[1] = scale;
}

void ParticleShaderOriented::setq(double *q)
{
	sides = (int)q[0];
	scale = q[1];
}

void ParticleShaderOriented::qname(char **qn)
{
	qn[0] = "Number of sides";
	qn[1] = "Size of disk";
}

void ParticleShaderOriented::drawShape(int i)
{
	Implicit *imp = 0;

	if (ARData)
		radius=ARData->r[i];

	gluDisk(quad,0,radius*scale,sides,sides);
	glDisable(GL_LIGHTING);
	glBegin(GL_TRIANGLES);
		glColor3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(radius, 0.0, 0.0);
		glVertex3f(0.0, 0.0, radius);

		glColor3f(1.0, 1.0, 1.0);
		glVertex3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, radius, 0.0);
		glVertex3f(0.0, 0.0, radius);
	glEnd();
	glEnable(GL_LIGHTING);
}

ParticleShaderOriented::~ParticleShaderOriented()
{
	gluDeleteQuadric(quad);
}
