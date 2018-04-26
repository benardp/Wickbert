/*
* @file SPHInitial.h
* @author Yuan Zhou
*/
#ifndef SPHINITIAL_H
#define SPHINITIAL_H

#include "Particles.h"
#include "ParticleBehavior.h"
#include "ParticleAcceleration.h"
#include "ParticleDensity.h"
#include "ParticleLocalityGrid.h"


class SPHInitial : public ParticleBehavior
{
public:
	ParticleDensity *density;
	ParticleAcceleration *acceleration;

	ParticleLocality* plocality;

	// create SPH attribute for Particles p, default constructor
	SPHInitial(Particles *ps = NULL);

	MAKE_PARTICLESTUFF_NAME();

	double mass;
	double density_kernel;
	double density_rest;

	double h2;
	double h9;
	double coeff;
	
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

	double W_poly6(const gmVector3 &r);

	virtual void attachAttributes();

	virtual void applyForce();

	virtual void cleanup();
};

#endif
