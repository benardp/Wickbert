/**
* Declaration of a silhouette detector behavior.
* @file SilhouetteDetector.h
* @date 19 May. 2005
* @author Matei N. Stroila (based on part of FeatureDetector by John C. Hart)
*/

#ifndef SILHOUETTEDETECTOR_H
#define SILHOUETTEDETECTOR_H


/** This class detects the silhouette (sampling with floaters) 
 * and creates particles to represent the silhouette 
 */
class SilhouetteDetector : public ParticleBehavior
{	
public:
	MAKE_PARTICLESTUFF_NAME();

	SilhouetteDetector(Particles *ps=NULL);
	~SilhouetteDetector();
	virtual void attachAttributes();
	// everything implemented in cleanup
	virtual void cleanup();
	void findSilhouettes();
protected:
	void particleAdded();
	void particleRemoved(unsigned int i);
	std::vector<gmVector3> lastN;
private:
	double anglevar;	///< Angle variance threshold
	std::string target;	///< name of Particles where new particles are created
	Particles *target_p;
	ParticlePosition *pos;
	ViewDependence *view;
	ImplicitInterrogator *impInt;
	ParticleNormal *orient;
	Implicit *imp;
};

#endif
