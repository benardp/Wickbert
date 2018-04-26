//RBF Propagation
//Mike Mullan, John C. Hart

#ifndef RBFPROPAGATION_H
#define RBFPROPAGATION_H

#include "ParticleBehavior.h"
#include "Particles.h"
#include "Surface/Implicit/Implicit.h"
#include "ImplicitInterrogator.h"
#include "ParticleOrientation.h"
#include "Surface/Implicit/Variational/RBFInterpolant.h"


//Moves a RBF surface by the rules of level set propagation
class RBFPropagation : public ParticleBehavior {

private:

	/// Interrogator for the surface to be deformed.
	ImplicitInterrogator* imp_int;

	/// Particle orientations.
	ParticleOrientation* p_orient;

	
public:

	MAKE_PARTICLESTUFF_NAME();

	/// Factor to scale the speed function
	double speedfac;
	double constfac;
	double gaussfac;
	double meanfac;

	/// speed function
	virtual double speed(gmVector3 point);

	int qlen() { return 4; }
	void getq(double *q) {q[0] = speedfac; q[1] = constfac; q[2] = gaussfac; q[3] = meanfac;}
	void setq(double *q) {speedfac = q[0]; constfac = q[1]; gaussfac = q[2]; meanfac = q[3];}
	void qname(char **qn) {qn[0] = "Speed factor"; qn[1] = "Constant rate";
						   qn[2] = "Gaussian curvature rate"; qn[3] = "Mean curvature rate";}

	

	void attachAttributes();

	/// Creates a surface adhesion attribute for Particles p.
	RBFPropagation(Particles *ps=NULL, const std::string& name=std::string("RBFPropagation"));
    
	/** Returns the particle orientation.
	*/
	ParticleOrientation* getOrientation() { return p_orient; }

	ImplicitInterrogator* getImplicitInterrogator() { return imp_int; }

	/// Apply the surface adhesion constraint.
	void applyConstraint();

	/// Update surface parameters based on actual dqdt values.
	void integrate();

	void propagateNormals(RBFInterpolant* rbf);  //Rotates the dipole normal vectors as the RBF surface is propagating
	double targetShape(gmVector3 point);
	gmVector3 computeGradF(gmVector3 point);  //The gradient of the speed function
	void adaptivity(RBFInterpolant* rbf);  //Adaptively add and delete RBF centers
	void RBFPropagation::tune(RBFInterpolant* rbf);


};

#endif