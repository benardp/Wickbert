/*
@file DuplicateParticlesAlongNormal.h
@author Elmar Eisemann
This attribute can be used to achieve RBF sampling representations

It initializes a particles object P based on a second particles object Q.
The entry are two scalars per particle o1,o2, the particle's position x and normal n.
P will contain twice the number of particles in Q. Namely x+o*n and x-o2*n.
Realize that if one stores a single offset it is possible to use it for both entries.
*/

#ifndef DUPLICATEPARTICLESALONGNORMAL_H
#define DUPLICATEPARTICLESALONGNORMAL_H

#include "../Particles.h"
#include "../ParticleAttribute.h"

class ParticlePosition;
class ParticleScalar;

class DuplicateParticlesAlongNormal : public ParticleAttribute {

	//unfortunately these two attributes do not inherit from ParticleVector.
	//I will change this later and then this attribute accepts ANY Vector attribute.
	ParticlePosition * _thisPosition;
	ParticleScalar * _thisScalar;

	
	ParticleAttribute * _positionAttribute;
	ParticleAttribute * _normalAttribute;
	ParticleAttribute * _scalarAttribute;

	//this offset is used, whenever no scalar per particle attribute has been provided.
	double _offset;

public:
	DuplicateParticlesAlongNormal(Particles *ps=NULL, const std::string& name=std::string("DuplicateParticleNormals"));
	void duplicateParticlesAlongNormal();
	
	MAKE_PARTICLESTUFF_NAME();
};

class ApplyDuplicationAlongNormalCallback : public PSParamButton::Callback
{
public:
	DuplicateParticlesAlongNormal *vp;
	ApplyDuplicationAlongNormalCallback (DuplicateParticlesAlongNormal *me) {vp = me;}
	virtual void onbuttonpress() {vp->duplicateParticlesAlongNormal();}
};

#endif
