/**
 * Declaration of SilhouetteFate.h
 * @file SilhouetteFate.h
 * @date February 16, 2002
 * @author Doug Nachand
 */

#ifndef SILHOUETTEFATE_H
#define SILHOUETTEFATE_H	

#include "Particles.h"
#include "ParticleBehavior.h"
#include "ParticleAge.h"
#include "ParticleRepulsion.h"
#include "ParticleBoundingBox.h"
#include "SilhouetteAdhesion.h"

/**
 * This is the class which governs particle birth and death as described by
 * Witkin and Heckbert with modifications to now consider the silhouette. 
 * Particles which are large and have low energy birth new particles and 
 * particles which are overcroweded are killed.
 */
class SilhouetteFate : public ParticleBehavior{

public:

	// attibutes
	AdaptiveRepulsionData* rep_data; ///< Witkin-Heckbert repulsion data.
	ParticleOrientation* p_orient;   ///< Particle orientations.
	ParticleAge* p_age;				///< Age of the silhouette particles
	SilhouetteAdhesion* sil_ad;	

	// behaviors
	ParticleBoundingBox* p_bounds;   ///< Bounding box.

	MAKE_PARTICLESTUFF_NAME();

	/// parameters
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

	double surfaceDiameter;

	double gamma; ///< Equilibrium speed (multiple of repulsion radius).
	double nu;    ///< Fraction of E_hat, for fissioning.
	double delta; ///< Fraction of sigma_hat, for death.
	int population;

	virtual void attachAttributes();

	SilhouetteFate(Particles *ps=NULL, ParticleBoundingBox *b=NULL,
		SilhouetteAdhesion *sa=NULL, const std::string& name=std::string("SilhouetteFate"));

	/// Allows for an awareness of surface diameter.
	double setSurfaceDiameter(double d);

	/*
	 *	Returns true if the number of particles on the silhouette
	 *	is at least fifty percent.
	 */
	bool fiftyPercent();

	void resetAge(){p_age->reset();}

	/// We determine particle fate during the cleanup step.
	void cleanup();

};

#endif

