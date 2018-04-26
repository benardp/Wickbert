/*
 *  ParabolicsAdhesion.cpp
 *  GPS
 *
 *  Created by Matei Stroila on 12/9/05.
 *  
 */

#include "ParabolicsAdhesion.h"
#include "Surface/Implicit/Operator/ZeroGaussianCurvature.h"

#include "ParticleNormal.h"
#include "ImplicitInterrogator.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"

REGISTER_PARTICLESTUFF(ParabolicsAdhesion,"Behavior:ParabolicsAdhesion");

ParabolicsAdhesion::ParabolicsAdhesion(Particles *ps)
: ParticleBehavior(ps, std::string("ParabolicsAdhesion"))
{
	new PSParamDouble(this,&phi,15.0,"k","penalty strength",
					  "Strength of the penalty constraint feedback term to "
					  "ensure the dynamic constraint keeps the particle at the zero level set");
	new Attached<ImplicitInterrogator>(this,&imp_int);
	new Attached<ParticleNormal>(this,&p_orient);
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
	ImplicitIsSet = false;
	imp1 = imp2 = NULL;
}

void ParabolicsAdhesion::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	ImplicitIsSet = false;
}

void ParabolicsAdhesion::applyConstraint()
{
	
	unsigned int i;			// index of current particle
	gmVector3 xi,				// position of particle i
		gradient,					// gradient of current particle
		secondNormal;					// gradient of the second surface
	double a,b,m, lambda, mu, Fi, Si,denom;
	
	if(!ImplicitIsSet) ImplicitIsSet = setImplicit();
	if(!ImplicitIsSet) return;
	
	// If we get to here, then we have an implicit surface to work with.
	
	DoubleArray dfdq(0.0,imp1->qlen());
	
	// Check size of desired parameter change array
	
	if (imp_int->dqdt.size() != dfdq.size())
	{
		imp_int->dqdt.resize(dfdq.size());
		imp_int->dqdt = 0.0;
	}
	
	for (i = 0; i < ps->size(); i++)
	{
		
		xi = position->getPosition(i);
		gradient = imp_int->grad(i);//imp1->grad(xi);
		
		p_orient->setNormal(i, gradient);
		p_orient->setGradMag2(i, gradient.lengthSquared());
		
		secondNormal = imp2->grad(xi);
		a = gradient.lengthSquared();
		b = secondNormal.lengthSquared();
		m = dot(gradient,secondNormal);
		denom = m*m - a*b;
		//if critical point do nothing
		if(gmIsZero(denom)) continue; 
		
		imp1->procq(xi, dfdq);
		Fi = phi* imp_int->proc(i)//imp1->proc(xi) 
			+ dot(gradient, velocity->v[i]);
		Si = phi* imp2->proc(xi) 
			+ dot(secondNormal, velocity->v[i]);
		lambda = -b*Fi + m*Si;
		mu = -a*Si + m*Fi;		
		velocity->v[i] += (-lambda * gradient - mu * secondNormal)/denom;		
	}
}

bool ParabolicsAdhesion::setImplicit(){
	
	if (!imp_int) return false;
	// get the implicit surface from the interrogator
	imp1 = imp_int->getImplicit();
	if (!imp1) return false;
	// the second implicit is the corresponding parabolic surface
	if(imp2) delete imp2;
	imp2 = new ZeroGaussianCurvature(imp1);
	return true;
}

/**
* Callback for particle addition.
 * @param i of the added particle.
 */
void ParabolicsAdhesion::particleAdded()
{
	if(!ImplicitIsSet) ImplicitIsSet = setImplicit();
	if(!ImplicitIsSet) return;
	int i = ps->size() - 1;
	//ELMAR DEBUG - how are we gonna do this?
	p_orient->setNormal(i, imp1->normal(position->getPosition(i)));
} 

ParabolicsAdhesion::~ParabolicsAdhesion()
{
	delete imp2;
}
