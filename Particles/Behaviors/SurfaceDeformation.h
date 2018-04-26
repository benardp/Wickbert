/**
 * Declaration of SurfaceDeformation.
 * @file SurfaceDeformation.h
 * @author John Hart, Ed Bachta
 */

#ifndef SURFACEDEFORMATION_H
#define SURFACEDEFORMATION_H

#include "ParticleBehavior.h"

#include "svd.h"
#include "tnt/tnt.h"

class ImplicitInterrogator;
class ParticlePosition;
class ParticleNormal;
class ParticleVelocity;
class ParticleScalar;
class Implicit;

/**
 * SurfaceDeformation allows the parameters of an implicit function
 * to be modified as a result of moving particles that are constrained
 * to lie on the surface. Our implementation follows the Witkin-Heckbert
 * model.
 *
 * SurfaceDeformation cannot be applied to a particle system that does
 * not yet have the ImplicitInterrogator attribute. This is to ensure
 * that it does not attempt to operate on a null surface pointer.
 */
class SurfaceDeformation : public ParticleBehavior {

private:

	SVD svd;

	/// Interrogator for the surface to be deformed.
	ImplicitInterrogator* imp_int;

	/// Particle orientations.
	ParticleNormal* p_orient;
	
	ParticlePosition *position;
	ParticleVelocity *velocity;

	/// Desired particle orientations.
	//ParticleDesiredOrientation* d_orient;

	void interpolate(Implicit *imp);

	bool solveCholesky(TNT::Matrix<double>& A,
		TNT::Vector<double>& x,
		TNT::Vector<double>& b);

public:

	MAKE_PARTICLESTUFF_NAME();

	/// Feedback coefficient.
	double phi;

	/** Desired function values at control points.
	 * Usually this is a vector of zeros,
	 * except for radial basis functions.
	 */
	ParticleScalar* scalar;
	std::vector<double> h;

	/** A vector indicating which procq elements are "flexible",
	 * which means they are affected by control particles.
	 *
	 * \todo flexible should be moved from SurfaceDeformation to Implicit.
	 */
	std::valarray<bool> flexible;

	void attachAttributes();

	/// Creates a surface adhesion attribute for Particles p.
	SurfaceDeformation(Particles *ps=NULL);

	void loadParticles();

	/// Apply the surface adhesion constraint.
	void applyConstraint();

	virtual void particleAdded();

	virtual void particleRemoved(unsigned int i);

};

#endif
