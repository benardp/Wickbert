// Commented this file out to remove depedency on ISM - colombi@uiuc.edu
#if 0

/*
* @file CriticalPointInterrogator.h
* @date October 8,2003
* @author Terry Fleury <tfleury@uiuc.edu>
* Transform the CriticalPoint code into a Particle Behavior
*/

#ifndef CRITICALPOINTINTERROGATOR_H
#define CRITICALPOINTINTERROGATOR_H

#include "Particles.h"
#include "ParticleBehavior.h"
#include "ParticleBoundingBox.h"
#include "ism/SearchCritical.h"

class ImplicitInterrogator;

class CriticalPointInterrogator : public ParticleBehavior
{
  private:
    /// Attribute for the Implicit Surface
    ImplicitInterrogator* impInt;

    /// Set to true to force recalculation of critcal points
    int calculateCriticalPoints;
    /// Pointer to the critical point particles
    Particles *target;
    /// The Newton solver to find 3D critical points
    SearchCritical criticalPointsFinder;
    /// Bounding box for the Newton solver
    ParticleBoundingBox* bounds;
    
  public:

    MAKE_PARTICLESTUFF_NAME();

    /// Default constructor
    CriticalPointInterrogator(Particles* ps=NULL,
      const std::string& name=std::string("CriticalPointInterrogator"));

    void attachAttributes();

    /// Parameter control for the user
    int qlen();
    void getq(double *q);
    void setq(double *q);
    void qname(char **qn);

    void createCriticalPointPS(Particles* ps);
    void cleanup();
    void findCriticalPoints();

};

#endif

#endif