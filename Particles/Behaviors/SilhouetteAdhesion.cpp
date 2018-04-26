/**
 * Implementation of the silhouette constraint.
 * @file SilhouetteAdhesion.h
 * @author Matei N. Stroila
 * @date 6/22/2005
 * @remarks 
 */
#include "SilhouetteAdhesion.h"

#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "ParticleNormal.h"
#include "ImplicitInterrogator.h"
#include "ViewDependence.h"
#include "ParticleAge.h"


REGISTER_PARTICLESTUFF(SilhouetteAdhesion,"Behavior:SilhouetteAdhesion");

SilhouetteAdhesion::SilhouetteAdhesion(Particles *ps)
: TwoImplicitIntersectionAdhesion(ps, std::string("SilhouetteAdhesion"))
{
	new PSParamDouble(this,&phi,15.0,"k","penalty strength",
					  "Strength of the penalty constraint feedback term to "
					  "ensure the dynamic constraint keeps the particle at the zero level set");
	new PSParamDouble(this,&pAgeThreshold, 1.0e-3,"pAgeThreshold","Age Threshold",
				"A particle won't age if the product of the values of the implicits defining the contour,"
				"evaluated at the particle position, is less than this value. ");
	new Attached<ImplicitInterrogator>(this,&imp_int);
	new Attached<ParticleNormal>(this,&p_orient);
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
	new Attached<ViewDependence>(this,&view);
}

void SilhouetteAdhesion::attachAttributes()
{
	ParticleBehavior::attachAttributes();
}
void SilhouetteAdhesion::initParameters()
{
	cameraPos = *(view->getCameraPosition());
}
void SilhouetteAdhesion::applyConstraint()
{
	
	initParameters();
	unsigned int i;			// index of current particle
	gmVector3 xi,				// position of particle i
		gradient,					// gradient of current particle
		silhNormal;					// gradient of the silhouette surface
	double a,b,m, lambda, mu, Fi, Si,denom;
	
	// get the implicit surface from the interrogator
	if (!imp_int) return;

	Implicit *imp = imp_int->getImplicit();
	
	// Implicit surface might not exist. In this case the interrogator exists but there is no Implicit for it to refer.
	if (!imp) return;
	
	// If we get to here, then we have an implicit surface to work with.
		
	for (i = 0; i < ps->size(); i++)
	{
		xi = position->getPosition(i);
		gradient = imp_int->grad(i);
		
		p_orient->setNormal(i, gradient);
		const double gradLengthSquared=gradient.lengthSquared();
		p_orient->setGradMag2(i, gradLengthSquared);
		
		// Particles with low gradients should stay put.
		if (gradLengthSquared < 0.001) continue;
		
		silhNormal = imp_int->hess(i)//getImplicit()->hess(xi) 
						* (cameraPos - xi) - gradient;
		a = gradLengthSquared;
		b = silhNormal.lengthSquared();
		m = dot(gradient,silhNormal);
		denom = m*m - a*b;
		//if critical point do nothing
		if(gmIsZero(denom)) continue; 
		
		Fi = phi*imp_int->proc(i)
			+ dot(gradient, velocity->v[i]);
		Si = dot(silhNormal, velocity->v[i])-phi*onSilhouetteSurfProc(xi,gradient);

		lambda = -b*Fi + m*Si;
		mu = -a*Si + m*Fi;		
		velocity->v[i] += (-lambda * gradient - mu * silhNormal)/denom;		

		pAge->resetParticleAge(i, fabs(imp_int->proc(i) * onSilhouetteSurfProc(xi,gradient)) < pAgeThreshold);
	}
}

double SilhouetteAdhesion::onSilhouetteSurfProc(const gmVector3& pos, const gmVector3& grad) const
{
		return dot(grad, pos-cameraPos );
}


double SilhouetteAdhesion::definingManifoldProc(const gmVector3 & position) const
{
	gmVector3 gradient=imp_int->getImplicit()->grad(position);
	return onSilhouetteSurfProc(position,gradient);
}

/**
* Callback for particle addition.
* @param i of the added particle.
*/
void SilhouetteAdhesion::particleAdded()
{
	Implicit *imp = imp_int->getImplicit();
	
	// Implicit interrogator may exist but might not refer to an Implicit
	if (!imp) return;
	
	int i = ps->size() - 1;
	//ELMAR DEBUG - how are we gonna do this?
	p_orient->setNormal(i, imp->normal(position->getPosition(i)));
} 