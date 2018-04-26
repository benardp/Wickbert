/**
 * Implementation of the shadow detector.
 * @file ShadowDetector.cpp
 * @date 27 May. 2005
 * @author Matei N. Stroila (based on part of FeatureDetector by John C. Hart)
 */
#include "ShadowDetector.h"

#include "ParticlePosition.h"
#include "ParticleNormal.h"
#include "LightPosition.h"
#include "Particles.h"

REGISTER_PARTICLESTUFF(ShadowDetector,"Behavior:ShadowDetector");

ShadowDetector::ShadowDetector(Particles *ps) 
	: ParticleBehavior(ps, std::string("ShadowDetector"))
{
	new PSParamDouble(this,&anglevar,9.0,
		"anglevar","angle variance",
		"Angle variance theshold to determine shadows.");

	new PSParamString(this,&target,"sil","sil","target particles",
		"Name of Particles collection within which to create new particles.");

	new Attached<ParticleNormal>(this,&orient);
	new Attached<ParticlePosition>(this,&pos);
	new Attached<LightPosition>(this,&lightPosAtt);
}

ShadowDetector::~ShadowDetector()
{
}

void ShadowDetector::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	if (!ps) return;
	target_p = ps->particleSystem->findParticles(target);
}

void ShadowDetector::cleanup()
{
	if (!target_p) return;

	lightPosAtt->getLightPosition(lightPos);

	gmVector3 lightPosition3(lightPos[0], lightPos[1], lightPos[2]);

	float dotmin = (float) cos((90.0 - anglevar)*M_PI/180.0);
	
	ParticlePosition *target_pos =
		target_p->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
	if (!target_pos) return;
	ParticleNormal *target_orient =							
		target_p->getAttribute<ParticleNormal>(std::string("ParticleNormal"));
	if (!target_orient) return;

	target_p->removeAll();
	gmVector3 ni, viewLighti;
	for (unsigned int i = 0; i < ps->size(); i++) {
		ni = orient->getNormal(i);
		viewLighti = lightPosition3 -  pos->getPosition(i);
		viewLighti.normalize();
		float doti = (float) dot(ni,viewLighti);
		if (fabs(doti) < fabs(dotmin)) {
			int added = target_p->addParticle();
			target_pos->setPosition(added,pos->getPosition(i));
			target_orient->setNormal(added,ni);
		}
	}
}
