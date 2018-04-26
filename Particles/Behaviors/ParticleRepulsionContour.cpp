/**
 * Implementation of the repulsion force for particles adhering to a contour.
 * @file ParticleRepulsionContour.cpp
 * @date July 17, 2006
 * @author Matei Stroila
 */


#include "ParticleRepulsionContour.h"

#include "pstools.h"
#include <math.h>
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "AdaptiveRepulsionData.h"
#include "Contours.h"


REGISTER_PARTICLESTUFF(ParticleRepulsionContour,"Behavior:ParticleRepulsionContour");

/**
* Creates a particle repulsion attribute for Particles p.
 * @param ps The owning particle system.
 */
ParticleRepulsionContour::ParticleRepulsionContour(Particles *ps) 
: ParticleBehavior(ps, std::string("ParticleRepulsionContour"))
{
	new PSParamDouble(this,&beta,10.0,
					  "beta","denominator fudge",
					  "Term added to denominator to prevent division by zero.");
	
	new PSParamDouble(this,&rho,15.0,
					  "rho","feedback","Strength of penalty force to keep particle on surface.");
	
	new Attached<AdaptiveRepulsionData>(this,&rep_data,"AdaptiveRepulsionData",
										"repdata","Repulsion Data",
										"Source of the parameters of the Gaussian repulsion model.");
	
	
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
	new Attached<Contours>(this,&contours);
}

void ParticleRepulsionContour::integrate()
{
	if (ps && rep_data)
		rep_data->integrate(ps->dt);
}

/** Need to update the grid in case
* particles have moved to different cells
*/
void ParticleRepulsionContour::cleanup()
{
}

/**
* Applies a force that repels nearby particles.
 * \todo Need to write this so the query particle need not be in the
 * same particles as its "neighbors" organized in ParticleLocality.
 */
void ParticleRepulsionContour::applyForce() {
		
	unsigned int i, j;
	double	rij2,	// distance squared between particle and one of its neighbors
		sig2,	// particle's repulsion radius squared
		Eij,	// Energy term between particle and one of its neighbors
		dDidr;	// Derivative of D[i] wrt radius
	gmVector3 rij;	// Vector between particles i and j.
		
	for (i = 0; i < ps->size(); ++i) {
		
		/* Find this particle's neighbors.
		*/
		
		Neighbors nb;
		nb.n1 = nb.n2 = i; 
		contours->getNeighbors(i,nb);
		unsigned int neighbors[2] = {nb.n1, nb.n2};
					
		rep_data->D[i] = 0.0;
		dDidr = 0.0;
		unsigned int *nbr;
		nbr = neighbors;
		// Apply the force based on all neighbors
		
		for (int k = 0; k < 2; k++) {
			j = *nbr;
			
			if (i == j) continue;	// Don't repel self
			
			rij = position->getPosition(i) - position->getPosition(j);
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
			
			/* Add energy into particle i's velocity (first half of WH (9))
				*/
			velocity->v[i] += rij*Eij;
			
			/* Add energy into particle j's velocity (second half of (9))
				* Negate because rij is in the opposite direction.
				* Note that the subtraction of the energy terms should actually
				* be an addition in (9). See Heckbert's erratum.
				*/
						
			velocity->v[j] -= (rep_data->r[j]*rep_data->r[j]/sig2)*rij*Eij;
			
			
			/* D[i] is the sum of the total energy.
				* Also defined as the portion of particle i's
				* energy directly affected by a change in its
				* own repulsion radius. [WH]
				*/
			rep_data->D[i] += Eij;
			
			/* dDidr is the derivative of Di wrt the repulsion radius.
				* See WH (13).
				*/
			dDidr += rij2 * Eij / (sig2*rep_data->r[i]);
			nbr++;
		} // end cache loop
		
		/* DiDot is the change in Di over time, used
			* as a feedback to control Di. See (10) in WH.
			*/
		double Didot = -rho*(rep_data->D[i] - rep_data->Ehat);
		
		/* dr[i] is the change in the particle i's radius.
			* Called \dot{sigma}^i in WH, see (12).
			*/
		rep_data->dr[i] = Didot / (dDidr + beta);		
	}
}


