/**
* @file SurfaceAdhesion.cpp
* Implementation of the surface constraint.
*/

#include "SurfaceAdhesion.h"

#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "ParticleNormal.h"
#include "ImplicitInterrogator.h"
#include "ParticleAge.h"

REGISTER_PARTICLESTUFF(SurfaceAdhesion,"Behavior:SurfaceAdhesion");

SurfaceAdhesion::SurfaceAdhesion(Particles *ps)
: ParticleBehavior(ps, std::string("SurfaceAdhesion"))
{
	new PSParamDouble(this,&phi,15.0,"k","penalty strength",
		"Strength of the penalty constraint feedback term to "
		"ensure the dynamic constraint keeps the particle at the zero level set");
	new PSParamDouble(this,&singTh,0.01,"singTh","singularity threshold",
					  "This number is the minimum allowed gradient length squared");
	new PSParamDouble(this,&pAgeThreshold, 1.0e-3,"pAgeThreshold","Age Threshold",
			"A particle won't age if the values of the implicit defining the surface,"
			"evaluated at the particle position, is less than this value. ");
	new Attached<ImplicitInterrogator>(this,&imp_int);
	new Attached<ParticleNormal>(this,&p_orient);
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
	new Attached<ParticleAge>(this,&pAge);
}

void SurfaceAdhesion::applyConstraint()
{
	unsigned int i;			// index of current particle
	gmVector3 xi,				// position of particle i
	gradient;					// gradient of current particle
	double val; 				// magnitude of constraint force
	std::vector<gmVector3> v;	// temporary velocity holder
	
	// This initializes the valarray to be filled with 0.0
	//DoubleArray dfdq(0.0,imp_int->getImplicit()->qlen());
	/** \bug SurfaceAdhesion::applyConstraint wants to access an implicit interrogator
	 * implicit surface, but the imp_int->getImplicit() pointer is bad when an implicit
	 * has been deleted. Probably a problem with using integer indices to connect to
	 * implicits.
	 */

	if (!imp_int) return;

	// get the implicit surface from the interrogator
	Implicit *imp = imp_int->getImplicit();

	// Implicit surface might not exist. In this case the interrogator exists but there is no Implicit for it to refer.
	if (!imp) return;

	// If we get to here, then we have an implicit surface to work with.

	DoubleArray dfdq(0.0,imp->qlen());

	// Check size of desired parameter change array
	if (imp_int->dqdt.size() != dfdq.size())
	{
		imp_int->dqdt.resize(dfdq.size());
		imp_int->dqdt = 0.0;
	}

	for (i = 0; i < ps->size(); i++)
	{
		
		xi = position->getPosition(i);
		gradient = imp_int->grad(i);//getImplicit()->grad(xi);

		const double gradLengthSquared=gradient.lengthSquared();

		//ELMAR DEBUG - Is this really a good place?
		p_orient->setNormal(i, gradient);
		p_orient->setGradMag2(i, gradient.lengthSquared());

#if 1
		// Particles with low gradients should stay put.
		if (gradLengthSquared < singTh)
			continue;
#endif

		imp_int->getImplicit()->procq(xi, dfdq);
		
		// Evaluate WH (5), portion of velocity in gradient direction
		val = dot(gradient, velocity->v[i]);
		
		// Compute procq . qdot which accounts for
		// motion of particle due to changing surface parameters
		if (dfdq.size() > 0 && imp_int->dqdt.size() > 0)
			val += (dfdq*imp_int->dqdt).sum();

		// Constrain will keep the function value constant
		// but due to numerical error, this constant can
		// drift away from zero. The phi parameter is the
		// spring constant of a good old penalty constraint
		// to bring the particle back to the surface if it
		// drifted.
		double proc = imp_int->proc(i);//getImplicit()->proc(xi);

	
		val += phi*proc;

		val /= gradLengthSquared;
		
		velocity->v[i] += -val*gradient;

		pAge->resetParticleAge(i,(fabs(proc) < pAgeThreshold));
		
	}	
}


/**
* Callback for particle addition.
* @param i of the added particle.
*/
void SurfaceAdhesion::particleAdded()
{
	Implicit *imp = imp_int->getImplicit();

	// Implicit interrogator may exist but might not refer to an Implicit
	if (!imp) return;

	//as usual, is this really a good place? - Elmar
	int i = ps->size() - 1;
	p_orient->setNormal(i, imp->normal(position->getPosition(i)));
} 
