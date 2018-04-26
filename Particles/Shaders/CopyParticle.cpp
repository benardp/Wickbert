#include "CopyParticle.h"
#include "ParticleSystem.h"

REGISTER_PARTICLESTUFF(CopyParticle,"Shader:CopyParticle");

CopyParticle::CopyParticle(Particles *ps)
	: ParticleShader(ps,std::string("CopyParticle"))
{
	new PSParamInt(this,&eventcode,-1,"eventcode","On Event",
		"Event codes: click = 1, doubleclick = 2, drag = 3, rightclick = 4...");
	new PSParamString(this,&targetparticles,"<empty>",
		"targetparticles","Target Particles",
		"Name of Particles collection into which to copy.");

	//new Attached<ParticlePosition>(this,&pos);
	new Attached<ParticlePosition>(this,&pos,"ParticlePosition","srcPos","Source Position","Where the particle positions are copied FROM.");
	new PSParamString(this, &dstPosName, "ParticlePosition","destPos","Destination Position","Where the particle positions are copied TO.");
}

void
CopyParticle::event(int e)
{
	if (e == eventcode) {
		Particles *tp = ps->particleSystem->findParticles(targetparticles);
		if (tp) {
			int i = ps->selectedParticle;
			int j = tp->addParticle();

			/* only copies position for now */

			//std::string("ParticlePosition")
			ParticlePosition *tpos= tp->getAttribute<ParticlePosition>(dstPosName);
			if(tpos)
				tpos->setPosition(j,pos->getPosition(i));
		}
	}
}
