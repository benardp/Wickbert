/**
 * Implementation of the silhouette detector.
 * @file SilhouetteDetector.h
 * @date 19 May. 2005
 * @author Matei N. Stroila (based on part of FeatureDetector by John C. Hart)
 */
#include "Particles.h"
#include "ParticleBehavior.h"
#include "ParticlePosition.h"
#include "ParticleNormal.h"
#include "ParticleSystem.h"
#include "libgm/gmVec3.h"
#include "ImplicitInterrogator.h"
#include "ViewDependence.h"
#include "SilhouetteDetector.h"

#ifndef M_PI
#define M_PI 3.1415
#endif


REGISTER_PARTICLESTUFF(SilhouetteDetector,"Behavior:SilhouetteDetector");

SilhouetteDetector::SilhouetteDetector(Particles *ps) 
	: ParticleBehavior(ps, std::string("SilhouetteDetector"))
{
	new PSParamDouble(this,&anglevar,9.0,
		"anglevar","angle variance",
		"Angle variance theshold to determine silhouette.");

	new PSParamString(this,&target,"sil","sil","target particles",
		"Name of Particles collection within which to create new particles.");
	new Attached<ParticlePosition>(this,&pos);
	new Attached<ViewDependence>(this,&view);
	new Attached<ImplicitInterrogator>(this,&impInt);
	new Attached<ParticleNormal>(this,&orient);
}

void SilhouetteDetector::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	if (!ps) return;
	target_p = ps->particleSystem->findParticles(target);
	
	lastN.resize(ps->size());
	for (unsigned int i = 0; i < ps->size(); i++)
		lastN[i] = orient->getNormal(i);
}

void SilhouetteDetector::particleAdded()
{
	lastN.push_back(gmVector3());
}

void SilhouetteDetector::particleRemoved(unsigned int i)
{
	lastN[i] = lastN.back();
	lastN.pop_back();
}

void SilhouetteDetector::cleanup()
{
	
	findSilhouettes();
	
}

void SilhouetteDetector::findSilhouettes()
{
	if (!target_p) return;
	float dotmin = (float)(cos((90.0 - anglevar)*M_PI/180.0));
	float dotminSing = (float)(cos(anglevar*M_PI/180.0));

	ParticlePosition *target_pos =
		target_p->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
	if (!target_pos) return;
	ParticleNormal *target_orient =	target_p->getAttribute<ParticleNormal>(std::string("ParticleNormal"));
	if (!target_orient) return;

	imp = impInt->getImplicit();
	ImplicitInterrogator *target_impInt =
		target_p->getAttribute<ImplicitInterrogator>(std::string("ImplicitInterrogator"));
	if (!target_impInt) return;
	target_impInt->setImplicit(imp);

	target_p->removeAll();
	gmVector3 viewVectori, ni, posi;
	
	for (unsigned int i = 0; i < ps->size(); i++) {
		posi = pos->getPosition(i);
		ni = orient->getNormal(i);
		viewVectori = (*view->getCameraPosition()) -  posi;
		viewVectori.normalize();
		float doti = (float)(dot(ni,viewVectori));
		float dotiSing = (float)(dot(ni,lastN[i]));
		if ((fabs(doti) < fabs(dotmin)) && ( dotiSing > dotminSing)) {
			int added = target_p->addParticle();
			target_pos->setPosition(added,posi);
			target_orient->setNormal(added,ni);
		}
		lastN[i] = ni;
	}//end for
}//end findSilhouettes()

SilhouetteDetector::~SilhouetteDetector()
{
}
