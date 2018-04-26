/**
* Implementation of the shadow constraint.
 * @file ShadowAdhesion.cpp
 * @date June 06, 2005
 * @author Matei Stroila
 */
#include "ShadowAdhesion.h"

#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "ParticleNormal.h"
#include "ImplicitInterrogator.h"
#include "LightPosition.h"
#include "ParticleAge.h"

REGISTER_PARTICLESTUFF(ShadowAdhesion,"Behavior:ShadowAdhesion");

ShadowAdhesion::ShadowAdhesion(Particles *ps)
: TwoImplicitIntersectionAdhesion(ps, std::string("ShadowAdhesion"))
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
	new Attached<LightPosition>(this,&light);
}

void ShadowAdhesion::attachAttributes()
{
	ParticleBehavior::attachAttributes();
}


void ShadowAdhesion::initParameters()
{
	gmVector4 lightPos4;
	light->getLightPosition(lightPos4);
	
	directionalLight=(lightPos4[3]==0.0);

	if(directionalLight)
		for(unsigned int j = 0; j <3 ; j++) lightPos3[j] = lightPos4[j];
	else
		for(unsigned int j = 0; j <3 ; j++) lightPos3[j] = lightPos4[j]/lightPos4[3];
}


void ShadowAdhesion::applyConstraint()
{

	initParameters();

	unsigned int i;			// index of current particle
	gmVector3 xi,				// position of particle i
		gradient,					// gradient of current particle
		shadowNormal;					// gradient of the shadow surface
	double a,b,m, lambda, mu, Fi, Si,denom;
	
	// get the implicit surface from the interrogator
	if (!imp_int) return;
	Implicit *imp = imp_int->getImplicit();
	
	// Implicit surface might not exist. In this case the interrogator exists but there is no Implicit for it to refer.
	if (!imp) return;
			
	for (i = 0; i < ps->size(); i++)
	{	
		xi = position->getPosition(i);
		gradient = imp_int->grad(i);//getImplicit()->grad(xi);
		
		p_orient->setNormal(i, gradient);
		p_orient->setGradMag2(i, gradient.lengthSquared());
		
		// Particles with low gradients should stay put.
		if (gradient.lengthSquared() < 0.001) continue;
		
		if(directionalLight)
			shadowNormal = imp_int->hess(i)//getImplicit()->hess(xi) 
							* (lightPos3);
		else
			shadowNormal = imp_int->hess(i)//getImplicit()->hess(xi) 
			* (lightPos3 - xi) - gradient;
		
		a = gradient.lengthSquared();
		b = shadowNormal.lengthSquared();
		m = dot(gradient,shadowNormal);
		denom = m*m - a*b;
		if(gmIsZero(denom)) continue; 
		
		Fi = phi*imp_int->proc(i)//imp_int->getImplicit()->proc(xi) 
				+ dot(gradient, velocity->v[i]);

		Si = phi*onShadowSurfProc(xi,gradient) + dot(shadowNormal, velocity->v[i]);
		
		lambda = -b*Fi + m*Si;
		mu = -a*Si + m*Fi;		
		velocity->v[i] += (-lambda * gradient - mu * shadowNormal)/denom;	
		pAge->resetParticleAge(i, fabs(imp_int->proc(i) * onShadowSurfProc(xi,gradient)) < pAgeThreshold);
	}
}

double ShadowAdhesion::definingManifoldProc(const gmVector3& pos) const
{
	const gmVector3 grad=imp_int->getImplicit()->grad(pos);
	return onShadowSurfProc(pos,grad);
}

double ShadowAdhesion::onShadowSurfProc(const gmVector3& pos, const gmVector3& grad) const
{
	if (directionalLight)
		return dot(grad, lightPos3);
	else
		return dot(grad, lightPos3-pos);
}

/**
* Callback for particle addition.
 * @param i of the added particle.
 */
void ShadowAdhesion::particleAdded()
{
	Implicit *imp = imp_int->getImplicit();
	
	// Implicit interrogator may exist but might not refer to an Implicit
	if (!imp) return;
	
	int i = ps->size() - 1;
	//ELMAR DEBUG - how are we gonna do this?
	//behaviors ask for attributes, but cannot ask for other behaviors.
	//thus who assures that the values in the attributes are what we would expect?
	p_orient->setNormal(i, imp->normal(position->getPosition(i)));
} 