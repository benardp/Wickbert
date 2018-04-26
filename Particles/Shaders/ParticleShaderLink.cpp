#include "ParticleShaderLink.h"
#include "ParticleLocality.h"
#include "AdaptiveRepulsionData.h"
#include "ParticlePosition.h"

REGISTER_PARTICLESTUFF(ParticleShaderLink,"Shader:ParticleShaderLink");

ParticleShaderLink::ParticleShaderLink(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderLink"))
{
	//new Attached<ParticleLocality>(this,&p_locality);
	new Attached<AdaptiveRepulsionData>(this,&ardata);
	new Attached<ParticleLocalityGrid>(this,&p_locality);
}

void ParticleShaderLink::drawPre()
{
	p_locality->update();
}

void ParticleShaderLink::drawParticle(int i)
{
	if (ps->selectedParticle != -1 && ps->selectedParticle != i) return;

	std::list<unsigned int> neighbors;
	std::list<unsigned int>::iterator nbr;
	int j;
		
	neighbors.clear();
	double queryRadius = ardata->sdmul * ardata->r[i];
	p_locality->getNeighbors(i,queryRadius,neighbors);

	gmVector3 vi = p_locality->position->getPosition(i);

	glBegin(GL_LINES);			
	for (nbr = neighbors.begin(); nbr != neighbors.end(); nbr++) {
		j = *nbr;
		gmVector3 vj = p_locality->position->getPosition(j);
		glVertex3d(vi[0],vi[1],vi[2]);
		glVertex3d(vj[0],vj[1],vj[2]);
	}
	glEnd();
}
