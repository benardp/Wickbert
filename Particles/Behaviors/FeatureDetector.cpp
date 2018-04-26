/**
 * Implementation of the feature detector.
 * @file FeatureDetector.cpp
 * @date 4 Feb. 2004
 * @author John C. Hart (based on SingularityAdhesion by Wen Su)
 */
#include "FeatureDetector.h"
#include "ParticlePosition.h"
#include "ParticleNormal.h"
#include "Particles.h"
#include "ParticleSystem.h"
#include "libgm/gmVec3.h"

#ifndef M_PI
#define M_PI 3.1415
#endif


REGISTER_PARTICLESTUFF(FeatureDetector,"Behavior:FeatureDetector");

FeatureDetector::FeatureDetector(Particles *ps) 
	: ParticleBehavior(ps, std::string("FeatureDetector"))
{
	new PSParamDouble(this,&anglevar,30.0,
		"anglevar","angle variance",
		"Angle variance theshold to determine feature.");
	new PSParamInt(this, &flip_number, 5, "flip number", "flip number", "The minimum number of flips for feature particle");

	new PSParamString(this,&target,"target","target","target particles",
		"Name of Particles collection within which to create new particles.");

	new Attached<ParticleNormal>(this,&orient);
	new Attached<ParticlePosition>(this,&pos);
}

void
FeatureDetector::attachAttributes()
{
	ParticleBehavior::attachAttributes();

	if (!ps) return;

	target_p = ps->particleSystem->findParticles(target);

	lastN.resize(ps->size());
	flips.resize(ps->size());

	for (size_t i = 0; i < ps->size(); i++)
		lastN[i] = orient->getNormal(i);
	for (size_t i = 0; i < ps->size(); i++)
		flips[i] = 0;
}

void FeatureDetector::particleAdded()
{
	lastN.push_back(gmVector3());
	flips.push_back(0);
}

void FeatureDetector::particleRemoved(unsigned int i)
{
	lastN[i] = lastN.back();
	lastN.pop_back();
	flips[i] = flips.back();
	flips.pop_back();
}

void
FeatureDetector::cleanup()
{
	if (!target_p) return;
	float dotmin = (float) cos(anglevar*M_PI/180.0);
	ParticlePosition *target_pos =
		target_p->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
	if (!target_pos) return;

	for (size_t i = 0; i < ps->size(); i++) {
		gmVector3 ni = orient->getNormal(i);
		float doti = (float) dot(ni,lastN[i]);
		if (doti < dotmin) {
			if( flips[i]++ == flip_number ) {
			int added = target_p->addParticle();
			target_pos->setPosition(added,pos->getPosition(i));
			}
		}
		lastN[i] = ni;
	}
}
