// Commented this file out to remove depedency on ISM - colombi@uiuc.edu
#if 0

/**
 * @file CriticalPointInterrogator.cpp
 * @date October 8, 2003
 * @author Terry Fleury <tfleury@uiuc.edu>
 */

#include "CriticalPointInterrogator.h"
#include "ImplicitInterrogator.h"
#include "ParticleSystem.h"
#include "ParticleShaderSphere.h"
#include "ParticlePosition.h"
#include "ParticleEigen.h"
#include "ParticleMaterial.h"

REGISTER_PARTICLESTUFF(CriticalPointInterrogator,"Behavior:CriticalPointInterrogator");

/// Default constructor.
CriticalPointInterrogator::CriticalPointInterrogator(Particles *ps, 
                                                     const std::string& name):
  ParticleBehavior(ps,name)
{
  impInt = NULL;
  target = NULL;
  calculateCriticalPoints = 0;

  // imp may be NULL at this point so move it to later
  criticalPointsFinder = SearchCritical(impInt->getImplicit());
}

void CriticalPointInterrogator::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	attachAttribute(impInt,std::string("ImplicitInterrogator"));
}

int CriticalPointInterrogator::qlen()
{
  return 1;
}

void CriticalPointInterrogator::getq(double *q)
{
  q[0] = calculateCriticalPoints;
}

void CriticalPointInterrogator::setq(double *q)
{
  calculateCriticalPoints = (int)q[0];
}

void CriticalPointInterrogator::qname(char **qn)
{
  qn[0] = "set 1 to recalculate critical points";
}

void CriticalPointInterrogator::createCriticalPointPS(Particles* ps)
{
  // If there is already a target or ps is NULL then return
  if (ps == NULL || target != NULL)
    return;

  target = new Particles(std::string("CriticalPoint"));

  if (!target) // Check for error creating the target 
    return;

  new ParticleEigen(target);
  new ParticlePosition(target);
	new ParticleMaterial(target);
  new ParticleShaderSphere(target);
  new ImplicitInterrogator(target, impInt->getImplicit());

  // Record it in the system
  if (ps->particleSystem)
    ps->particleSystem->particles.push_back(target);
}

void CriticalPointInterrogator::cleanup()
{
  // Do nothing unless the user asked to calculate critical points
  if (!calculateCriticalPoints)
    return;
  calculateCriticalPoints = 0;

  createCriticalPointPS(ps);  // Should only run once

  // createCriticalPointsPS should create target, but maybe an error?
  if (target == NULL)
    return;
  else
    target->clear();

  if (impInt && impInt->getImplicit())
	 criticalPointsFinder = SearchCritical(impInt->getImplicit());
  else
	 return;
  // Do the actual work of finding the critical points with the Newton solver
  findCriticalPoints();

  // Once we have found the critical points, copy them into a particle system
  for (unsigned int idx = 0; idx < criticalPointsFinder.fRoots.size(); ++idx)
    {
      target->addParticle();
      int lastIdx = target->size() - 1;

      ParticlePosition *targetPosition = 
        target->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
      ParticleEigen *targetEigen = 
        target->getAttribute<ParticleEigen>(std::string("ParticleEigen"));
			ParticleMaterial *targetMaterial =
				target->getAttribute<ParticleMaterial>(std::string("ParticleMaterial"));

      targetPosition->setPosition(lastIdx,
          gmVector3(criticalPointsFinder.fRoots[idx][0],
                    criticalPointsFinder.fRoots[idx][1],
                    criticalPointsFinder.fRoots[idx][2]));
      targetEigen->setEigenVectors(lastIdx,criticalPointsFinder.fVectors[idx]);
      targetEigen->setEigenValues(lastIdx,criticalPointsFinder.fValues[idx]);
      targetEigen->setPositiveEigenValues(lastIdx,
        (int)criticalPointsFinder.fKinds[idx]);

			switch((int)criticalPointsFinder.fKinds[idx])
				{
					case 0: targetMaterial->setDiffuse(
										lastIdx,gmVector4(0.2,1.0,0.2,1.0)); break; // MAX
					case 1: targetMaterial->setDiffuse(
										lastIdx,gmVector4(0.2,1.0,1.0,1.0)); break; // SADDLE2
					case 2: targetMaterial->setDiffuse(
										lastIdx,gmVector4(1.0,1.0,0.2,1.0)); break; // SADDLE2
					case 3: targetMaterial->setDiffuse(
										lastIdx,gmVector4(1.0,0.2,0.2,1.0)); break; // MIN
				}
    }
}

void CriticalPointInterrogator::findCriticalPoints()
{
	// this should really use attachAttribute.
  //bounds = ps->getBehavior<ParticleBoundingBox>(std::string("ParticleBoundingBox"));
  //if (bounds == NULL)
  //  bounds = new ParticleBoundingBox(ps);

  //criticalPointsFinder.setBoxBounds(bounds->getBounds());
  //criticalPointsFinder.fUseNewton = false;
  //criticalPointsFinder.search();
  //criticalPointsFinder.classify();
}


#endif