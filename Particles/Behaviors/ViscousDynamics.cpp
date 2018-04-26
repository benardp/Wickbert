/**
 * @file ViscousDynamics.cpp (originally VelocityIntegrator.cpp) 
 * @author Wen Su
 */

#include "ViscousDynamics.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "ParticleScalar.h"

REGISTER_PARTICLESTUFF(ViscousDynamics,"Behavior:ViscousDynamics");

ViscousDynamics::ViscousDynamics(Particles *ps, const std::string& name) 
	: ParticleBehavior(ps, name)
{
	new PSParamDouble(this,&dt,0.03,
		"dt","Time Step","Time step used for Euler integration.");

	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
//	new Attached<ParticleScalar>(this,&speed);
	new Attached<ParticleVelocity>(this,&history);
}

void ViscousDynamics::integrate()
{
	for(unsigned int i=0;i<ps->size();i++)
	{
		position->setPosition(i,position->getPosition(i)+dt*velocity->v[i]);
//		speed->setScalar(i, (velocity->v[i]).length());
		history->v[i] = dt * velocity->v[i];
		velocity->v[i] = gmVector3();
	}
}
