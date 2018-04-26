/**
 * @file ParticleShader.cpp
 * @author Wen Su
 */

#include "ParticleShader.h"
#include "ParticleMesh.h"
// #include "ParticleOrientation.h"
// #include "ParticleDesiredOrientation.h"
#include "ParticleMaterial.h"
//#include "ParticleVisibility.h"
#include "AdaptiveRepulsionData.h"

ParticleShader::ParticleShader(Particles *ps, const std::string& name)
	: ParticleStuff(ps,name)
{

	iterations = 0;
	pre_time = draw_time = post_time = 0;

}

void ParticleShader::setParticleSystem(Particles *new_ps)
{
	if (!new_ps)
		return;
	ps = new_ps;
	addMeTo(ps);
}

void ParticleShader::attachAttributes()
{
	attachedattributes.attach();
}

void ParticleShader::removeMe()
{
	if (ps==NULL)
		return;

	bool move = false;
	for(unsigned int i = 0; i < ps->shaders.size(); i++)
	{
		if (move)
			ps->shaders[i-1] = ps->shaders[i];
		if(ps->shaders[i]->getName() == name)
			move = true;
	}
	if (move)
		ps->shaders.resize(ps->shaders.size()-1);
}

bool ParticleShader::moveAfter(ParticleStuff *stuff)
{
	unsigned int i,j,k;

	ParticleShader *shad = dynamic_cast<ParticleShader *>(stuff);
	if (!shad) return false;
	for (i = 0; i < ps->shaders.size(); i++)
		if (ps->shaders[i] == this) break;
	if (i == ps->shaders.size()) return false; /* and confused because I should be in the shaders! */
	for (j = 0; j < ps->shaders.size(); j++)
		if (ps->shaders[j] == shad) break;
	if (j == ps->shaders.size()) return false; /* stuff came from another particles's shaders */

	/* now need to move element i to element after j */

	if (i < j) {
		for (k = i; k < j; k++)
			ps->shaders[k] = ps->shaders[k+1];
		ps->shaders[j] = this;
	} else {
		for (k = i; k > j+1; k--)
			ps->shaders[k] = ps->shaders[k-1];
		ps->shaders[j+1] = this;
	}
	return true;
}

void ParticleShader::addMeTo(Particles *new_ps)
{
	// So help me, if new_ps is null, I'm bringing the app down.
	assert(new_ps);

	// allow shaders with the same name
	new_ps->shaders.push_back(this);

	attachAttributes();
}

