/*
Declaration of SingularityInterrogator
@file SingularityInterrogator.h
@author Wen Su
This take the initial guess of SingularityAdhesion, refine it by minimizing the square norm of the gradient.
*/

#ifndef SINGULARITYINTERROGATOR_H
#define SINGULARITYINTERROGATOR_H

#include "Particles.h"
#include "ParticleBehavior.h"

class ImplicitInterrogator;

class SingularityInterrogator : public ParticleBehavior
{
private:
	// parameter
	// scale the gradient
	double gradientConst;
	double feedbackConst;
	// interval grad size
	double initIntervalSize;
	double maxLength;

public:
	
	MAKE_PARTICLESTUFF_NAME();
	
	/// parameters to control for the user
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

	/// default constructor
	SingularityInterrogator(Particles* ps=NULL, const std::string& name = std::string("SingularityInterrogator"));
	
	ImplicitInterrogator* impInt;

	// move opposite of gradient
	void applyForce();

	void attachAttributes();

};

#endif

