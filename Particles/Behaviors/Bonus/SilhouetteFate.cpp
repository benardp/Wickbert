/**
 * Implementation of ParticleFate.
 * @file SilhouetteFate.cpp 
 * @date November 28, 2001
 * @author Ed Bachta
 */

#include "SilhouetteFate.h"
#include "pstools.h"

// Witkin-Heckbert default values
#define WH_GAMMA 4.0
#define WH_NU 0.2
#define WH_DELTA 0.7

REGISTER_PARTICLESTUFF(SilhouetteFate,"Behavior:SilhouetteFate");

/**
 * Creates a SilhouetteFate object. First it looks for the required
 * attributes in the particle system, and then it sets default 
 * constants.
 */
SilhouetteFate::SilhouetteFate(Particles *ps, ParticleBoundingBox *b, SilhouetteAdhesion *sa,
	const std::string& name) : ParticleBehavior(ps, std::string("SilhouetteFate"))
{
	gamma = WH_GAMMA;
	nu = WH_NU;
	delta = WH_DELTA;

	surfaceDiameter = 1;
	population = 100;

	// these two are behaviors
	p_bounds=b;
	sil_ad=sa;
}

/**
 *
 * @todo This behavior uses data from the SilhouetteAdhesion behavior.
 * The shared data should be made into an Attribute.
 */
void SilhouetteFate::attachAttributes()
{
	ParticleBehavior::attachAttributes();

	attachAttribute(rep_data,std::string("AdaptiveRepulsionData"));
	attachAttribute(p_age,std::string("ParticleAge"));
	attachAttribute(p_orient,std::string("ParticleOrientation"));
	attachAttribute(p_bounds,std::string("ParticleBoundingBox"));

	p_age->reset();
	setSurfaceDiameter(surfaceDiameter);

	// get behaviors
	sil_ad=ps->getBehavior<SilhouetteAdhesion>(std::string("SilhouetteAdhesion"));

}

/** 
 * Allows for an awareness of surface rep_data->diameter.
 * @param d New surface rep_data->diameter.
 * @returns New desired readius.
 */
double SilhouetteFate::setSurfaceDiameter(double d)
{
	if (d > 0.0) 
	{
		rep_data->diameter = d;
		rep_data->sigma_hat = rep_data->diameter/4.0;
		rep_data->sigma_max = rep_data->diameter/2.0;
	} 
	return rep_data->sigma_hat;
}

/**
 * We determine particle fate during the cleanup step.
 */
void SilhouetteFate::cleanup()
{
	setSurfaceDiameter(surfaceDiameter);

	unsigned int i,j;
	gmVector3 N;
	//if (sil_ad->cameraMoved)
	//{
	//	p_age->reset();
	//	sil_ad->cameraMoved = false;
	//}
	
	// bool splitOnlyOnSilhouette = fiftyPercent();

	// Iterate backward to preserve indexing while removing particles
	// ps->size is unsigned int, so 0-1 can be bad
	for (j=ps->size(); j>0; j--) 
	{
		i=j-1;
		N = p_orient->getNormal(i);
		N.normalize();
		/**
		 * Determine whether or not particle is to die:
		 * 1) Particle velocity is small compared to nominal repulsion radius
		 * 2) Repulsion radius is smaller than minimum repulsion radius
		 * 3) Biased random test succeds
		 */
		if (((velocity->v[i].length() < gamma*rep_data->r[i]) &&
			(rep_data->r[i] < delta*rep_data->sigma_hat) &&
			(psRandom() > rep_data->r[i] / (delta*rep_data->sigma_hat))))
		{
			// kill this particle
			ps->removeParticle(i);
		}
 	}

	// Iterate forward to new size and split particles
	for (i=0; i<ps->size(); i++) 
	{
		N = p_orient->getNormal(i);
		N.normalize();
		/**
		* Determine whether or not to split:
		* (Condition 1:  particle is near equilibrium) AND
		* (Condition 2a: repulsion radius is huge OR
		*  Condition 2b: particle is adequately energized and repulsion radius
		*                is higher than nominal)
		*/
		if ( ((velocity->v[i].length() < gamma*rep_data->r[i]) &&
			((rep_data->r[i] > rep_data->sigma_max))))
			//|| 
			//((rep_data->D[i] > nu*rep_data->Ehat) && 
			//(rep_data->r[i] > rep_data->sigma_hat)))) )
		{
      
			// Caclulate perturbation
			double pert = 0.5*rep_data->r[i];
      
			// Generate tangent vectors orthogonal to the normal
			gmVector3 tgt1 = randOrtho(p_orient->getNormal(i));
			tgt1.normalize();

			gmVector3 tgt2 = cross(p_orient->getNormal(i), tgt1);
			tgt2.normalize();

			gmVector3 disp;

			int z=0;
			// Create a random displacement in bounds
			do 
			{
				disp = tgt1 * (2.0 * pert * psRandom() - pert) +
					tgt2 * (2.0 * pert * psRandom() - pert);
					pert *= 1.01;
				z++;

			} while ((p_bounds) && (!p_bounds->inBounds(position->getPosition(i) + disp)) && z < 50);

			// The new particle is added
			unsigned int j=ps->size();
			if (j<population)
			{
				ps->addParticle();
				position->setPosition(j,position->getPosition(i) + disp);
			}
			// Reduce the repulsion radii for both particles
			rep_data->r[i] /= sqrt(2.0);
			rep_data->r[j] = rep_data->r[i];
			p_age->t[j] = p_age->t[i];
		}
		else
			p_age->t[i] += ps->dt;
	}
}

int SilhouetteFate::qlen()
{
	return 2;
}
void SilhouetteFate::getq(double *q)
{
	q[0] = surfaceDiameter;
	q[1] = population;
}
void SilhouetteFate::setq(double *q)
{
	surfaceDiameter = q[0];
	population = q[1];
}
void SilhouetteFate::qname(char **qn)
{
	qn[0] = "Surface Diameter";
	qn[1] = "Population Limit";
}
