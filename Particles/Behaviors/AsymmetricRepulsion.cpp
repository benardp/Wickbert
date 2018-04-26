#include "AsymmetricRepulsion.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "pstools.h"


REGISTER_PARTICLESTUFF(AsymmetricRepulsion, "Behavior:AsymmetricRepulsion");

AsymmetricRepulsion::AsymmetricRepulsion(Particles* ps) : ParticleBehavior(ps, std::string("AsymmetricRepulsion"))
{
	new PSParamDouble(this, &factor, 1.0, "factor", "factor", "factor");
	new Attached<AdaptiveRepulsionData>(this,&rep_data,"AdaptiveRepulsionData",
		"repdata","Repulsion Data",
		"Source of the parameters of the Gaussian repulsion model.");

	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
	new PSParamString(this, &source, "source", "source", "source particles",
		"Name of Particles collection from which to apply the asymmetric repulsions.");
}

void AsymmetricRepulsion::applyForce()
{
	//	calculate repulsion among own particles
	//	This calculation calculates the self repulsion among own particles,
	//	it also calculate the differential of repulsion radius, but
	//	these differentials are not integrated until the integrate step,
	//	so the module below is still using the same repulsion radius
	//	as in the self-repulsion calculation.
	//	Another problem is that the source particles must have
	//	the AdaptiveRepulsionData attribute, but if it doesn't, 
	//	we are using own particles' radii as a estimation,
	//	afterall, the source particles are not supposed to move, and the values
	//	of their repulsion radii are just a scale factor of their significance.

	if( !source_p ) return;

	//	get the positions of the source particles
	ParticlePosition* s_position;
	source_p->findAttribute(s_position);

	AdaptiveRepulsionData* s_repdata;
	source_p->findAttribute(s_repdata);

	double rij2, sigi2, sigj2, Eij, Eji;
	gmVector3 rij;
	for(unsigned int i = 0; i < ps->size(); i++) {
		for(unsigned int j = 0; j < source_p->size(); j++) {
			rij = position->getPosition(i) - s_position->getPosition(j);
			rij2 = rij.lengthSquared();
			sigi2 = rep_data->r[i] * rep_data->r[i];
			if( s_repdata )
				sigj2 = s_repdata->r[j] * s_repdata->r[j];
			else
				sigj2 = sigi2;

			Eij = rep_data->alpha * fastExp(-0.5 * rij2 / sigi2);
			Eji = rep_data->alpha * fastExp(-0.5 * rij2 / sigj2);



//			Eij *= 1.0 - rij2 / (sig2 * rep_data->sdmul * rep_data ->sdmul);
			//	velocity constributions.
			velocity->v[i] += factor * rij * Eij;// + (rij * Eij * sigi2 / sigj2);
		}
	}
}

void AsymmetricRepulsion::attachAttributes()
{
	ParticleBehavior::attachAttributes();

	if( !ps ) return;
	source_p = ps->particleSystem->findParticles(source);
}
