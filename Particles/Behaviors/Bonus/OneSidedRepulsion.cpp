/**
 * @file OneSidedRepulsion.cpp
 * @date Jan 23, 2006
 * @author John C. Hart
 */

#include "ParticleLocalityGrid.h"
#include "ParticleRepulsion.h"
#include "pstools.h"
#include <math.h>
#include "ParticlePosition.h"
#include "ParticleVelocity.h"

REGISTER_PARTICLESTUFF(OneSidedRepulsion,"Behavior:OneSidedRepulsion");

/**
 * Creates a particle repulsion attribute for Particles p.
 * @param ps The owning particle system.
 */
ParticleRepulsion::ParticleRepulsion(Particles *ps) 
	: ParticleBehavior(ps, std::string("ParticleRepulsion"))
{
	new PSParamDouble(this,&beta,10.0,
		"beta","denominator fudge",
		"Term added to denominator to prevent division by zero.");

	new PSParamDouble(this,&rho,15.0,
		"rho","feedback","Strength of penalty force to keep particle on surface.");

	new Attached<AdaptiveRepulsionData>(this,&rep_data,"AdaptiveRepulsionData",
		"repdata","Repulsion Data",
		"Source of the parameters of the Gaussian repulsion model.");

	new Attached<ParticleLocality>(this,&otherlocality,"ParticleLocality",
		"locality","Particle Locality",
		"Indicates for each particle which particles are repelled.");

	new Attached<ParticlePosition>(this,&position,"ParticlePosition","pos","Repelling Positions","Positions of particles causing the motion.");
	new Attached<ParticlePosition>(this,&otherposition,"ParticlePosition","otherpos","Repelled Positions","Positions of particles being moved.");
	new Attached<ParticleVelocity>(this,&othervelocity,"ParticleVelocity","othervel","Repelled Velocities","Velocities of particles being moved.");
}

void ParticleRepulsion::integrate()
{
}

/** Need to update the grid in case
 * particles have moved to different cells
 */
void ParticleRepulsion::cleanup()
{
}
	
/**
 * Applies a force that repels nearby particles.
 * \todo Need to write this so the query particle need not be in the
 * same particles as its "neighbors" organized in ParticleLocality.
 */
void ParticleRepulsion::applyForce() {

	otherlocality->update();

	unsigned int i,k;
	int j;
	double	rij2,	// distance squared between particle and one of its neighbors
			sig2,	// particle's repulsion radius squared
			Eij,	// Energy term between particle and one of its neighbors
			dDidr;	// Derivative of D[i] wrt radius
	gmVector3 rij;	// Vector between particles i and j.
	
	std::list<unsigned int> neighbors;
	std::list<unsigned int>::iterator nbr;
		
	for (i = 0; i < ps->size(); ++i) {
		
		/* Find this particle's neighbors.
		 * neighbors is a vector of the indices
		 * of neighboring particles.
		 * Get the indices of neighbors that
		 * are as much as three standard deviations
		 * away.
		 * Note that the actual result of getNeighbors depends
		 * on the specifics of the locality object.
		 */

		neighbors.clear();
		double queryRadius = (rep_data->sdmul*rep_data->r[i]);
		p_locality->getNeighbors(i, queryRadius, neighbors);
		
		rep_data->D[i] = 0.0;
		dDidr = 0.0;

		// Apply the force based on all neighbors
		for (nbr = neighbors.begin(); nbr != neighbors.end(); nbr++) {
			j = *nbr;

			if (i == j) continue;	// Don't repel self

			rij = position->getPosition(i) - otherposition->getPosition(j);
			rij2 = rij.lengthSquared();
			sig2 = rep_data->r[i] * rep_data->r[i];
			
			/* Compute energy contribution
			 * This is the adaptive repulsion formula in WH 4.3
			 */
			Eij = rep_data->alpha*fastExp(-0.5*rij2/sig2);

			/* Scale Eij by distance so it goes to zero when |rij|
			 * nears sdmul*r[i].
			 * This keeps the particles stable when they are near
			 * the boundary of influence of other particles.
			 */
			Eij *= 1.0 - rij2/(sig2*rep_data->sdmul*rep_data->sdmul);
			
			/* Add energy into particle j's velocity (second half of (9))
			 * Negate because rij is in the opposite direction.
			 * Note that the subtraction of the energy terms should actually
			 * be an addition in (9). See Heckbert's erratum.
			 */

			othervelocity->v[j] -= (rep_data->r[j]*rep_data->r[j]/sig2)*rij*Eij;
		}
	}
}


