/*
@file ParticleBehavior
@author John Hart, Ed Bachta, Wen Su
Particles is too big, it separates into ParticleAttribute and ParticleBehavior
*/

#include "Particles.h"
#include "ParticleStuff.h"
#include "ParticleBehavior.h"

//#include "ParticlePosition.h"
//#include "ParticleMesh.h"
//#include "ParticleVelocity.h"

// Default contructor
ParticleBehavior::ParticleBehavior(Particles* ps, const std::string &name)
:ParticleStuff(ps,name)
{
	iterations = 0;
	total_time = 0;
}

void ParticleBehavior::addMeTo(Particles *new_ps)
{
	// So help me, if new_ps is null, I'm bringing the app down.
	assert(new_ps);

	int i = new_ps->findBehavior(name);
	if (i == -1)
		new_ps->behaviors.push_back(this);
	else
		new_ps->behaviors[i] = this;

	// attach all attributes needed
	attachAttributes();
}

void ParticleBehavior::attachAttributes()
{
#if 0
	// if ParticleMesh exists use it.
	ParticleMesh *pm=ps->getAttribute<ParticleMesh>(std::string("ParticleMesh"));
	if (pm)
		position=pm;
	else
		attachAttribute(position,std::string("ParticlePosition"));

	attachAttribute(velocity,std::string("ParticleVelocity"));
#endif

	attachedattributes.attach();
}

void ParticleBehavior::removeMe()
{
	if (ps==NULL)
		return;

	bool move = false;
	for(unsigned int i = 0; i < ps->behaviors.size(); i++)
	{
		if (move)
			ps->behaviors[i-1] = ps->behaviors[i];
		if(ps->behaviors[i]->getName() == name)
			move = true;
	}
	if (move)
		ps->behaviors.resize(ps->behaviors.size()-1);
}

bool ParticleBehavior::moveAfter(ParticleStuff *stuff)
{
	unsigned int i,j,k;

	ParticleBehavior *beha = dynamic_cast<ParticleBehavior *>(stuff);
	if (!beha) return false;
	for (i = 0; i < ps->behaviors.size(); i++)
		if (ps->behaviors[i] == this) break;
	if (i == ps->behaviors.size()) return false; /* and confused because I should be in the behaviors! */
	for (j = 0; j < ps->behaviors.size(); j++)
		if (ps->behaviors[j] == beha) break;
	if (j == ps->behaviors.size()) return false; /* stuff came from another particles's behaviors */

	/* now need to move element i to element after j */

	if (i < j) {
		for (k = i; k < j; k++)
			ps->behaviors[k] = ps->behaviors[k+1];
		ps->behaviors[j] = this;
	} else {
		for (k = i; k > j+1; k--)
			ps->behaviors[k] = ps->behaviors[k-1];
		ps->behaviors[j+1] = this;
	}
	return true;
}

