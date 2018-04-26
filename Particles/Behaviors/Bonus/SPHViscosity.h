/*
* @file SPHViscosity.h
* @author Yuan Zhou
*/
#ifndef SPHVISCOSITY_H
#define SPHVISCOSITY_H

#include "Particles.h"
#include "ParticleBehavior.h"
#include "ParticleAcceleration.h"
#include "ParticleVelocity.h"
#include "ParticleDensity.h"
#include "ParticleLocalityGrid.h"


class SPHViscosity : public ParticleBehavior
{
public:
	ParticleDensity *density;
	ParticleAcceleration *acceleration;

	ParticleLocality* plocality;

	// create SPH attribute for Particles p, default constructor
	SPHViscosity(Particles *ps = NULL);

	MAKE_PARTICLESTUFF_NAME();

	double mass;
	double kernel_h;
	double viscosity;

	double h2;
	double h6;
	double coeff;
	
	double ddW_visco(gmVector3 &r);

	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

	virtual void attachAttributes();

	virtual void applyForce();

	virtual void cleanup();
};

#endif
