/**
* Implementation of a curve adhesion behavior
 * @file CurveAdhesion.cpp
 * @date 9/27/2005
 * @author Matei N. Stroila
 * @remarks
 */

#include "CurveAdhesion.h"
#include "ParticleNormal.h"
#include "ImplicitInterrogator.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"


REGISTER_PARTICLESTUFF(CurveAdhesion,"Behavior:CurveAdhesion");

CurveAdhesion::CurveAdhesion(Particles *ps)
: ParticleBehavior(ps, std::string("CurveAdhesion"))
{
	new PSParamDouble(this,&phi,15.0,"k","penalty strength",
					  "Strength of the penalty constraint feedback term to "
					  "ensure the dynamic constraint keeps the particle at the zero level set");
	new Attached<ImplicitInterrogator>(this,&imp_int1);
	new Attached<ImplicitInterrogator>(this,&imp_int2);
	new Attached<ParticleNormal>(this,&p_orient);
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
}

void CurveAdhesion::attachAttributes()
{
	ParticleBehavior::attachAttributes();
}

void CurveAdhesion::applyConstraint()
{
	
	unsigned int i;			// index of current particle
	gmVector3 xi,				// position of particle i
		gradient,					// gradient of current particle
		secondNormal;					// gradient of the second surface
	double a,b,m, lambda, mu, Fi, Si,denom;
	DoubleArray dfdq(0.0,imp_int1->getImplicit()->qlen());
	
	// Check size of desired parameter change array
	
	if (imp_int1->dqdt.size() != dfdq.size())
	{
		imp_int1->dqdt.resize(dfdq.size());
		imp_int1->dqdt = 0.0;
	}
	
	for (i = 0; i < ps->size(); i++)
	{
		
		gradient = imp_int1->grad(i);
		
		p_orient->setNormal(i, gradient);
		p_orient->setGradMag2(i, gradient.lengthSquared());
		
		// Particles with low gradients should stay put.
		if (gradient.lengthSquared() < 0.001) continue;
		
		secondNormal = imp_int2->grad(i);
		a = gradient.lengthSquared();
		b = secondNormal.lengthSquared();
		m = dot(gradient,secondNormal);
		denom = m*m - a*b;
		//if critical point do nothing
		if(gmIsZero(denom)) continue; 
		

		//I dont think we want to store parameter derivatives in ImplicitInterrogator
		//but this could change one day - Elmar
		xi = position->getPosition(i);
		imp_int1->getImplicit()->procq(xi, dfdq);
		Fi = phi*imp_int1->proc(i) + dot(gradient, velocity->v[i]);
		Si = phi*imp_int2->proc(i) + dot(secondNormal, velocity->v[i]);
		lambda = -b*Fi + m*Si;
		mu = -a*Si + m*Fi;		
		velocity->v[i] += (-lambda * gradient - mu * secondNormal)/denom;		
	}
}


/**
* Callback for particle addition.
 * @param i of the added particle.
 */
void CurveAdhesion::particleAdded()
{
	int i = ps->size() - 1;
	p_orient->setNormal(i, imp_int1->normal(i));
} 


