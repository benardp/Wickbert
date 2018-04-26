/*
* @file SPHPressure.h
* @author Yuan Zhou
*/
#ifndef SPHPressure_H
#define SPHPressure_H

#include "Particles.h"
#include "ParticleBehavior.h"
#include "ParticleAcceleration.h"
#include "ParticleDensity.h"
#include "ParticleLocalityGrid.h"

class SPHPressure : public ParticleBehavior
{
public:
	ParticleDensity *density;
	ParticleAcceleration *acceleration;

	ParticleLocality* plocality;

	// create SPH attribute for Particles p, default constructor
	SPHPressure(Particles *ps = NULL);

	MAKE_PARTICLESTUFF_NAME();

	double mass;
	double gasC;
	double pressure_kernel;
	double density_rest;

	double h2;
	double h6;
	double coeff;
	
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

	gmVector3 dW_spiky(gmVector3 &r);
	
	virtual void attachAttributes();

	virtual void applyForce();

	virtual void cleanup();
};

#endif
