/**
 * Implementation of the specular highlights adhesion behavior 
 * @file SpecularAdhesion.cpp
 * @date 06/26/2005
 * @author Matei N. Stroila
 * @remarks It works only for orthogonal view and directional light.
 */
#include "SpecularAdhesion.h"

#include "ImplicitInterrogator.h"
#include "LightPosition.h"
#include "ViewDependence.h"
#include "ParticleNormal.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "ParticleAge.h"


REGISTER_PARTICLESTUFF(SpecularAdhesion,"Behavior:SpecularAdhesion");

SpecularAdhesion::SpecularAdhesion(Particles *ps)
: TwoImplicitIntersectionAdhesion(ps, std::string("SpecularAdhesion"))
{
	new PSParamDouble(this,&phi,15.0,"k","penalty strength",
					  "Strength of the penalty constraint feedback term to "
					  "ensure the dynamic constraint keeps the particle at the zero level set");
	new PSParamDouble(this,&shine,0.9,"shine","size of specular",
					  "The arccos of this number determines the specular shadow curve size");
	new PSParamDouble(this,&pAgeThreshold, 1.0e-3,"pAgeThreshold","Age Threshold",
					  "A particle won't age if the product of the values of the implicits defining the contour,"
					  "evaluated at the particle position, is less than this value. ");
	new Attached<ImplicitInterrogator>(this,&imp_int);
	new Attached<ParticleNormal>(this,&p_orient);
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
	new Attached<LightPosition>(this,&light);
	new Attached<ViewDependence>(this,&view);
}

void SpecularAdhesion::attachAttributes()
{
	ParticleBehavior::attachAttributes();
}

void SpecularAdhesion::initParameters()
{
	gmVector3 usedCameraPos = *(view->getCameraPosition());
	gmVector4 lightPos4;
	light->getLightPosition(lightPos4);
	gmVector3 lightPos3(lightPos4[0],lightPos4[1],lightPos4[2]);
	halfVector = (usedCameraPos.normalize() + lightPos3.normalize());
	halfVectorLength =  halfVector.length(); 
}

void SpecularAdhesion::applyConstraint()
{
	initParameters();
	
	unsigned int i;			// index of current particle
	gmVector3 xi,				// position of particle i
		gradienti,					// gradient of current particle
		specSurfGradienti;					// gradient of the specular surface
	double a,b,m, lambda, mu, Fi, Si,denom;
	
	// get the implicit surface from the interrogator
	if (!imp_int) return;
	Implicit *imp = imp_int->getImplicit();
	
	// Implicit surface might not exist. In this case the interrogator exists but there is no Implicit for it to refer.
	if (!imp) return;


	//make sure we don't divide by 0
	if(gmIsZero(halfVectorLength)) return;

	for (i = 0; i < ps->size(); i++)
	{		
		xi = position->getPosition(i);
		gradienti = imp_int->grad(i);//imp->grad(i)
		double gradiLengthSquare=gradienti.lengthSquared();
		double gradiLength = sqrt(gradiLengthSquare);
		
		//if the gradient is almost 0 at the position of this particle do not update the particle velocity 
		if(gmIsZero(gradiLength)) continue;
				
	
		p_orient->setNormal(i, gradienti);
		p_orient->setGradMag2(i, gradiLengthSquare);
		
		gmMatrix3 Hessi = imp_int->hess(i);//imp->hess(xi);
		
		specSurfGradienti = (Hessi * halfVector) -  shine *(halfVectorLength * (1/gradiLength) * (Hessi * gradienti));	
		 
		a = gradienti.lengthSquared();
		b = specSurfGradienti.lengthSquared();
		m = dot(gradienti,specSurfGradienti);
		denom = m*m - a*b;
		if(gmIsZero(denom)) continue; 
				
		Fi = phi*imp_int->proc(i)
		+ dot(gradienti, velocity->v[i]);
		Si = dot(specSurfGradienti, velocity->v[i]) - phi* onSpecSurfProc(gradienti,gradiLength);
			 			
		lambda = -b*Fi + m*Si;
		mu = -a*Si + m*Fi;		
		velocity->v[i] += (-lambda * gradienti - mu * specSurfGradienti)/denom;		

		pAge->resetParticleAge(i, fabs(imp_int->proc(i) * onSpecSurfProc(gradienti,gradiLength)) < pAgeThreshold);
	}
}

double SpecularAdhesion::onSpecSurfProc(const gmVector3 & gradient, double gradLength) const
{
	return (shine * gradLength *halfVectorLength  - dot(gradient, halfVector));
}
double SpecularAdhesion::definingManifoldProc(const gmVector3 & position) const
{
	gmVector3 gradient(imp_int->getImplicit()->grad(position));
	double gradLength=gradient.length();

	return onSpecSurfProc(gradient,gradLength);
}

/**
* Callback for particle addition.
 * @param i of the added particle.
 */
void SpecularAdhesion::particleAdded()
{
	Implicit *imp = imp_int->getImplicit();
	
	// Implicit interrogator may exist but might not refer to an Implicit
	if (!imp) return;
	
	int i = ps->size() - 1;
	//ELMAR DEBUG - same as in other classes, why is this class responsable for it...
	p_orient->setNormal(i, imp->normal(position->getPosition(i)));
} 
