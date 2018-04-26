/*
* @file SPHTension.h
* @author Yuan Zhou
*/
#ifndef SPHTENTION_H
#define SPHTENTION_H

#include "Particles.h"
#include "ParticleBehavior.h"
#include "ParticleAcceleration.h"
#include "ParticleDensity.h"
#include "ParticleLocalityGrid.h"


class SPHTension : public ParticleBehavior
{
public:
	ParticleDensity *density;
	ParticleAcceleration *acceleration;

	ParticleLocality* plocality;

	// create SPH attribute for Particles p, default constructor
	SPHTension(Particles *ps = NULL);

	MAKE_PARTICLESTUFF_NAME();

	double mass;
	double kernel_h;
	double tensionC;
	double normal_thresh;

	double h2;
	double h9;
	double coeff;

	gmVector3 dW_poly6(gmVector3 &r);
	double ddW_poly6(gmVector3 &r);

	
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

	virtual void attachAttributes();

	virtual void applyForce();

	virtual void cleanup();
};

#endif
