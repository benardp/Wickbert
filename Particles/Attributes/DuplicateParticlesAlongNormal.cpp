/**
 * @file ParticleVelocity.cpp
 * @author Elmar Eisemann
 */

#include "DuplicateParticlesAlongNormal.h"
#include "ParticlePosition.h"
#include "ParticleNormal.h"
#include "ParticleScalar.h"

REGISTER_PARTICLESTUFF(DuplicateParticlesAlongNormal,"Attribute:DuplicateParticlesAlongNormal");

/**
 * Add particle age to a system of particles.
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */
DuplicateParticlesAlongNormal::DuplicateParticlesAlongNormal(Particles *ps, const std::string& name)
	: ParticleAttribute(ps, name) 
{

	new Attached<ParticlePosition>(this, &_thisPosition);
	new Attached<ParticleScalar>(this, &_thisScalar);	

	new PSAttrRefParam(this, &_positionAttribute, "invalid", "position", "position for particles", "This value is used to define number and initial position of the particles (accepts position attribute)");

	new PSAttrRefParam(this, &_normalAttribute, "invalid", "normal", "normal for particles", "This value is used to define a normal at the initial positions (accepts normal attribute, attention need same number of particles!)");

	new PSParamDouble(this, &_offset, 0.01,"alt. offset", "alternative offset", "this offset is used globally if the per particle offset did not lead to an acceptable result");

	new PSAttrRefParam(this, &_scalarAttribute, "invalid", "offset", "offset for particles", "This value is used to define the distance travelled along and in the opposite direction of the normal (accepts any scalar attribute, attention need same number of particles!)");


	new PSParamButton(this, new ApplyDuplicationAlongNormalCallback(this),"duplicate","duplicate along normal",
		"We add particles that are positioned at .");
}

void DuplicateParticlesAlongNormal::duplicateParticlesAlongNormal()
{
	ParticlePosition * position=dynamic_cast<ParticlePosition*>(_positionAttribute);
	if (position==0)
		return;

	ParticleNormal * normal=dynamic_cast<ParticleNormal*>(_normalAttribute);
	if (normal==0)
		return;

	ParticleScalar * scalar=dynamic_cast<ParticleScalar*>(_scalarAttribute);

	unsigned int nbParticles = position->ps->size();
	if (normal->ps->size()!=nbParticles)
		return;


	if ((scalar==0)||(scalar->ps->size()!=nbParticles))
	{
		//duplicating with the fixed epsilon
		for (unsigned int i=0;i<nbParticles;++i)
		{
			unsigned int added=ps->addParticle();
			_thisPosition->setPosition(added, position->getPosition(i)+normal->getNormal(i)*_offset);
			_thisScalar->setScalar(added,_offset);
			added = ps->addParticle();
			_thisPosition->setPosition(added, position->getPosition(i)-normal->getNormal(i)*_offset);
			_thisScalar->setScalar(added,-_offset);
		}
	}
	else
	{
		for (unsigned int i=0;i<nbParticles;++i)
		{
			double offset = scalar->getScalar(i);
			unsigned int added=ps->addParticle();
			_thisPosition->setPosition(added, position->getPosition(i)+normal->getNormal(i)*offset);
			_thisScalar->setScalar(added, offset);
			added = ps->addParticle();
			_thisPosition->setPosition(added, position->getPosition(i)-normal->getNormal(i)*offset);
			_thisScalar->setScalar(added, -offset);
		}
	}
}
