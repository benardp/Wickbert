/**
* Declaration of a shadow detector behavior.
* @file ShadowDetector.h
* @date 27 May. 2005
* @author Matei N. Stroila (based on part of FeatureDetector by John C. Hart)
*/

#ifndef SHADOWDETECTOR_H
#define SHADOWDETECTOR_H

#include "ParticleBehavior.h"

class ParticleNormal;
class LightPosition;
class ParticlePosition;
class Particles;


/** This class detects the shadow of LIGHT0 (sampling with floaters) 
 * and creates particles to represent the shadow 
 */
class ShadowDetector : public ParticleBehavior
{	
public:
	MAKE_PARTICLESTUFF_NAME();

	ShadowDetector(Particles *ps=NULL);
	~ShadowDetector();
	virtual void attachAttributes();
	// everything implemented in cleanup
	virtual void cleanup();
private:
	double anglevar;	///< Angle variance threshold
	std::string target;	///< name of Particles where new particles are created
	Particles *target_p;
	ParticlePosition *pos;
	ParticleNormal *orient;
	float *lightPositionf;
	double *lightPositiond;
	LightPosition *lightPosAtt;
	gmVector4 lightPos;
};

#endif
