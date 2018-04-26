/**
 * Declaration of the repulsion force for particles adhering to a contour.
 * @file ParticleRepulsionContour.h
 * @date July 17, 2006
 * @author Matei Stroila
 */


#ifndef PARTICLEREPULSIONCONTOUR_H
#define PARTICLEREPULSIONCONTOUR_H

#include "ParticleBehavior.h"

class ParticlePosition;
class ParticleVelocity;
class AdaptiveRepulsionData;
class Contours;

/**
* This is the repulsion force described by Witkin and Heckbert S96
* which evenly distributes particles across the zero-set of an
* implicit function. The equations are derived assuming a single
* component of finite area.
*
* Holds two parameters. The parameter beta is a fudge added to the denominator
* in case it is zero. The parameter rho is a factor controlling the rate at
* which energy approaches a desired energy level.
*/
class ParticleRepulsionContour : public ParticleBehavior
{	
public:
	MAKE_PARTICLESTUFF_NAME();
	
	AdaptiveRepulsionData* rep_data;
	ParticlePosition* position;
	ParticleVelocity* velocity;
	Contours* contours;
	
	double beta; ///< Constant to prevent divide-by-zero.
	double rho;  ///< Feedback coefficient.
	
	/// Creates a particle repulsion attribute for Particles p.
	ParticleRepulsionContour(Particles *ps=NULL);
	
	/// Applies a force that repels nearby particles.
	void applyForce();
	
	// update repulsion data
	void integrate();
	
	// call update on locality
	void cleanup();	
};

#endif
