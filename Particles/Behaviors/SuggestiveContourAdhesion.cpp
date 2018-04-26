/**
* Implementation of the suggestive contour adhesion behavior
* @file SuggestiveContourAdhesion.cpp
* @date 23 June. 2005
* @author Matei N. Stroila
*/
#include "SuggestiveContourAdhesion.h"

#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "ParticleNormal.h"
#include "ImplicitInterrogator.h"
#include "ParticleCreation.h"
#include "ViewDependence.h"
#include "LightPosition.h"
#include "ParticleAge.h"


REGISTER_PARTICLESTUFF(SuggestiveContourAdhesion,"Behavior:SuggestiveContourAdhesion");

SuggestiveContourAdhesion::SuggestiveContourAdhesion(Particles *ps)
: TwoImplicitIntersectionAdhesion(ps, std::string("SuggestiveContourAdhesion"))
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
	
void SuggestiveContourAdhesion::initParameters()
{
	cameraPos = *(view->getCameraPosition());
}

void SuggestiveContourAdhesion::attachAttributes()
{
	ParticleBehavior::attachAttributes();
}

void SuggestiveContourAdhesion::applyConstraint()
{	
	initParameters();
	pCreation = ps->getBehavior<ParticleCreation>(std::string("ParticleCreation"));
	
	//start optimize particles' velocities
	unsigned int i;			// index of current particle
	gmVector3 xi,				// position of particle i
		gradient,					// gradient of the impilicit at the position of the current particle
		suggNormal;					// gradient of the suggestive surface
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
		gradient = imp_int->grad(i);//imp_int->getImplicit()->grad(xi);
		a = gradient.lengthSquared();
		
		p_orient->setNormal(i, gradient);
		p_orient->setGradMag2(i, a);
		
		// Particles with low gradients should stay put.
		if (a < 0.001)
			continue;
		gmVector3 v;
		v = cameraPos - xi; //view vector
	
		//check if the suggestive contour is defined
		if(!isDefined(v,p_orient->getNormal(i))) continue;
		
		//we have to think about this, if we precalculate all hessians it is probably slower than evaluating just some... - ELMAR
		gmMatrix3 Hess = imp_int->getImplicit()->hess(xi) ; //hessian
		gmVector3 Hessv = Hess * v;
		gmVector3 ugly(
					   dot((imp_int->getImplicit()->hessi(xi,0) * v),v),
					   dot((imp_int->getImplicit()->hessi(xi,1) * v),v),
					   dot((imp_int->getImplicit()->hessi(xi,2) * v),v)
					   );
		suggNormal = ugly - 2 * Hessv ;
		b = suggNormal.lengthSquared();
		// Particles with low sugg surfaces gradients should stay put.
		if (b < 0.001)
			continue;
		
		m = dot(gradient,suggNormal);
		denom = m*m - a*b;
		//if singularity, continue
		if(gmIsZero(denom)) continue;
	
		Fi = phi*imp_int->proc(i)
				+ dot(gradient, velocity->v[i]);
		Si = phi*dot(Hessv,v) + dot(suggNormal, velocity->v[i]);
		
		lambda = -b*Fi + m*Si;
		mu = -a*Si + m*Fi;		
		
		velocity->v[i] += (-lambda * gradient - mu * suggNormal)/denom;		

		pAge->resetParticleAge(i, fabs(imp_int->proc(i) * dot(Hessv,v)) < pAgeThreshold);
	}
}

double SuggestiveContourAdhesion::onSuggSurfProc(const gmVector3& pos, const gmVector3& Hessv) const
{
	return dot(Hessv, pos-cameraPos );
}

double SuggestiveContourAdhesion::definingManifoldProc(const gmVector3 & position) const
{
	gmVector3 Hessv = imp_int->getImplicit()->hess(position) * (position-cameraPos ) ;
	return onSuggSurfProc(position,Hessv);
}

bool SuggestiveContourAdhesion::isDefined(gmVector3 dir, gmVector3 normal) const
{
	bool defined = true;
	if(fabs(cross(dir.normalize(), normal).lengthSquared())<gmEPSILON) 
	{
		defined = false;
		pCreation->limitSize = 1;
	}
	return defined;
} 

/**
 * Callback for particle addition.
 * @param i the added particle.
 */
void SuggestiveContourAdhesion::particleAdded()
{
	Implicit *imp = imp_int->getImplicit();
	
	// Implicit interrogator may exist but might not refer to an Implicit
	if (!imp) return;
	
	int i = ps->size() - 1;
	p_orient->setNormal(i, imp->normal(position->getPosition(i)));
} 
