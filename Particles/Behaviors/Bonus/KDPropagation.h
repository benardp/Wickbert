
/**
 * Declaration of KDPropagation.
 * @file SurfaceDeformation.h
 * @author John Hart, Ed Bachta
 */

#ifndef KDPROPAGATION_H
#define KDPROPAGATION_H

//#include "Particles.h"
//#include "ParticleBehavior.h"
////#include "Implicit/Implicit.h"
//#include "ImplicitInterrogator.h"
//#include "ParticleOrientation.h"
//#include "svd.h"

#include "SurfacePropagation.h"
#include "KDTree/KDTree.h"

//class Implicit;
//class ParticlePosition;
//class ParticleVelocity;
//
//// TNT includes
//#include "tnt/tnt.h"
//#include "tnt/vec.h"
//#include "tnt/cholesky.h"
//#include "tnt/lu.h"
//#include "tnt/cmat.h"
//#include "tnt/transv.h"
//#include "tnt/trisolve.h"
#include <vector>
using namespace std;

struct kdpointt
{
	double pos[3];
	double norm[3];
};

class KDCompt
{
public:
	double operator() (const kdpointt& pt, int i) {
		return pt.pos[i];
	}
};

class KDTraitst
{
public:
	typedef kdpointt		point_type;
	typedef double		coord_type;
	typedef	KDCompt		get_component;
	static const int dim = 3;
};



/**
 * KDPropagation allows the parameters of an implicit function
 * to be modified as a result of a "level set" style evolution.
 * The particles are used to evaluate surface properties like curvature
 * to compute a per-particle speed function. We then solve a WH
 * control-particle system to change the parameters of the implicit
 * surface such that the particles move in the desired direction.
 *
 * Need to use per-particle error (difference in motion derived from change
 * in parameters and desired motion due to speed function) to cue the
 * subdivision of particles.
 */
class KDPropagation : public SurfacePropagation {

protected:

	///// Interrogator for the surface to be deformed.
	//ImplicitInterrogator* imp_int;

	//ParticlePosition *position;
	//ParticleVelocity *velocity;

	///// Singular Value Decomposition solver unit
	//SVD svd;

	///// Particle orientations.
	//ParticleOrientation* p_orient;

	///// Solves a system of equations using Cholesky factorization.
	//bool solveCholesky(TNT::Matrix<double>& A,
	//				 TNT::Vector<double>& x,
	//				 TNT::Vector<double>& b);
	//bool solveLU(TNT::Matrix<double>& A, TNT::Vector<double>& b);

	void load_point_cloud();

public:

	MAKE_PARTICLESTUFF_NAME();

	///// Factor to scale the speed function
	//double speedfac;
	//double constfac;
	//double gaussfac;
	//double meanfac;

	/// speed function
	virtual double speed(int i);

	/// Creates a surface adhesion attribute for Particles p.
	KDPropagation(Particles *ps=NULL, std::string name=std::string("KDPropagation"));

//	/** Returns the particle orientation.
//	*/
//	ParticleOrientation* getOrientation() { return p_orient; }
//
//	ImplicitInterrogator* getImplicitInterrogator() { return imp_int; }
//
//	/// Apply the surface adhesion constraint.
//	void applyConstraint();
//	void applyConstraintOLD();
//
//	/// Update surface parameters based on actual dqdt values.
//	void integrate();
//
//	/// Surface adhesion doesn't handle fissioning or death
////	void cleanup();

	xntools::ds::KDTree<KDTraitst> kd;
	vector<kdpointt> points;
	std::string ptcloud_file;
	bool loaded;
	bool reload;
};

#endif