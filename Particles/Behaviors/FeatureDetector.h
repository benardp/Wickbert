/**
* Declaration of a feature detector behavior.
* @file FeatureDetector.h
* @date 4 Feb. 2005
* @author John C. Hart (based on part of SingularityAdhesion by Wen Su)
*/

#ifndef FEATUREDETECTOR_H
#define FEATUREDETECTOR_H

#include "Particles.h"
#include "ParticleBehavior.h"

class ParticlePosition;
class ParticleNormal;

/** This class keeps a history of the particle orientations, and if there
 * is a lot of variance, it creates a new particle in another collection
 */
class FeatureDetector : public ParticleBehavior
{	
public:
	MAKE_PARTICLESTUFF_NAME();

	double anglevar;	///< Angle variance threshold

	std::string target;	///< name of Particles where new particles are created
	Particles *target_p;

	ParticlePosition *pos;
	ParticleNormal *orient;

	std::vector<gmVector3> lastN;
	int		flip_number;
	std::vector<int>	flips;

	/// Creates a particle repulsion attribute for Particles p.
	FeatureDetector(Particles *ps=NULL);

	virtual void attachAttributes();
	virtual void particleAdded();
	virtual void particleRemoved(unsigned int i);

	// everything implemented in cleanup
	virtual void cleanup();	
};

#endif
