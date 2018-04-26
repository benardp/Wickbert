#include "OrientParticle.h"
#include "ParticleNormal.h"
#include "ParticlePosition.h"

REGISTER_PARTICLESTUFF(OrientParticle,"Shader:OrientParticle");

OrientParticle::OrientParticle(Particles *ps)
	:ParticleShader(ps,std::string("OrientParticle"))
{
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleNormal>(this,&orientation);
}

void OrientParticle::drawParticle(int i)
{
	gmVector3 p = position->getPosition(i);
	glTranslatef((GLfloat) p[0], (GLfloat) p[1], (GLfloat) p[2]);

	gmMatrix4 rotMat = orientation->getMatrix(i);

	GLdouble mat[16];
	rotMat.copyTo(mat);

	glMultMatrixd(mat);	
}
