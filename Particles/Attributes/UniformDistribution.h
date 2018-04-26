/*
@file UniformDistribution.h
@author Elmar Eisemann
This attribute can be used to sample any kind of implicit function
The idea is that this attribute puts particles uniformly throughout the bounding box.
Now the choice of an appropriate shader gives the wanted effect

First I wanted to create a special shader to print out an implicit function, 
but this is more general. The shader could still be used on other particles that are not uniformly
distributed, for example to verify their values/normals.
And this uniform distribution could be used to initialize a particle system.
*/

#ifndef UNIFORMDISTRIBUTION_H
#define UNIFORMDISTRIBUTION_H

#include "../Particles.h"
#include "../ParticleAttribute.h"

class ParticlePosition;

class UniformDistribution : public ParticleAttribute {

	//unfortunately these two attributes do not inherit from ParticleVector.
	//I will change this later and then this attribute accepts ANY Vector attribute.
	ParticlePosition * _thisPosition;
	
	int _nbParticles[3];
	ParticleAttribute * _boundingBoxAttributeRef;
	bool _sampleMaxBorder;

public:
	UniformDistribution(Particles *ps=NULL, const std::string& name=std::string("UniformDistribution"));
	void redistributeUniformly();
	
	MAKE_PARTICLESTUFF_NAME();
};

class ApplyUniformRedistributionCallback : public PSParamButton::Callback
{
public:
	UniformDistribution *ud;
	ApplyUniformRedistributionCallback  (UniformDistribution *me) {ud = me;}
	virtual void onbuttonpress() {ud->redistributeUniformly();}
};

#endif
