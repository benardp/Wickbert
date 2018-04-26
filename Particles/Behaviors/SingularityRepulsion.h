/*
@file SingularityRepulsion.h
@author Wen Su
This class allows singularities to repel floaters in a different particle system.
*/

#ifndef SINGULARITYREPULSION_H
#define SINGULARITYREPULSION_H

#include "ParticleBehavior.h"

class ImplicitInterrogator;
class Particles;
class ParticlePosition;
class ParticleVelocity;
class ParticleRepulsion;

#define REPULSION_NEIGHBORS 6

class SingularityRepulsion : public ParticleBehavior
{	
protected:

	/// target system
	Particles *target_p;
	std::string target;
	/// target surface
	//do you need this ? - Elmar
	//ImplicitInterrogator *impInt;
	
	/// use the ParticleRepulsion system in that target
	ParticleRepulsion *repulsion;
	
	ParticlePosition *pos;
	ParticleVelocity *vel;
	
	double targetRepelConst;
	double selfRepelConst;
		
public:

	MAKE_PARTICLESTUFF_NAME();

	/// parameters
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

	/// default constructor
	SingularityRepulsion(Particles *ps=NULL/*, Particles *t=NULL, ParticleRepulsion *r=NULL*/);

	/// Applies a force that repels nearby particles.
	void applyForce();

	void repelTarget();
	void repelSelf();
	
	void attachAttributes();

};

#endif
