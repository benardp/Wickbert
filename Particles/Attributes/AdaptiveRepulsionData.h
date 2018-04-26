/** 
* Declaration of AdaptiveRepulsionData.
* @file AdaptiveRepulsionData.h
* @author Ed Bachta
*/

#ifndef __ADAPTIVEREPULSIONDATA_H__
#define __ADAPTIVEREPULSIONDATA_H__

#include "ParticleAttribute.h"

/**
* AdaptiveRepulsionData contains the parameters r, dr, and D along with
* constants which are needed to implement the Witkin-Heckbert model of
* adaptive particle repulsion. It is normally used by the ParticleRepulsion
* and ParticleFate behaviors to that effect.
*/
class AdaptiveRepulsionData : public ParticleAttribute {

public:

	MAKE_PARTICLESTUFF_NAME();

	// System attributes
	double sigma_hat; ///< Desired repulsion radius.
	double sigma_max; ///< Maximum repulsion radius.
	double diameter;	///< Surface diameter.
	double alpha; 	   ///< Repulsion amplitude.
	double Ehat;		 ///< Desired energy.

	/** # of std. dev.'s away that we expect repulsion to be zero
	 * Defaults to 3 standard deviations.
	 */
	double sdmul;

	// Per particle attributes

	std::vector<double> r;	///< Repulsion radii.
	std::vector<double> dr;	///< Change in repulsion radii.
	std::vector<double> D;	///< Particle energies.

	AdaptiveRepulsionData(Particles* ps=NULL,const std::string& name = std::string("AdaptiveRepulsionData"));

	virtual void setParticleSystem(Particles *);
	
	void clear();

	/** 
	* Integrate data over one time step.
	* @param dt Change in time to apply.
	*/
	void integrate(double dt);

	/// Callback for particle removal.
	virtual void particleRemoved(unsigned int i);
	
	/// Callback for particle addition.
	virtual void particleAdded();

	char *qtip(int i) { switch(i) {
		case 0:	/* sigma_hat */	return "Desired distance between particles.";
		case 1:	/* sigma_max */	return "Maximum distance between particles.";
		case 2: /* diameter */	return "Estimate of the diameter of the surface.";
		case 3:	/* alpha */		return "Amplitude of the repulsion force. Proportional to particle motion response time, "
									"but if too large may cause stiffness";
		case 4:	/* Ehat */		return "Desired energy of the system. Too much energy leads to jittery particles. "
									"Too little causes the particle system to appear sluggish.";
		case 5: /* sdmul */		return "Repulsion force is a Gaussian. The sdmul parameter is the number of "
									"standard deviations away from a particle that its repulsion force is "
									"tapered to zero.";
		default: return "";
	}}

	char *tip() { return "The AdaptiveRepulsionData particle attribute collects data shared "
						 "by several behaviors that collectively keep particles spaced "
						 "evenly.";
	}
};

#endif
