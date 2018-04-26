/**
Declaration of SingularityAdhesion
@file SingularityAdhesion.h
@author Wen Su
This class records the distance traveled for each particle and the angle of change for the normal to detect possible singularities
*/

#ifndef SINGULARITYADHESION_H
#define SINGULARITYADHESION_H

#include "Particles.h"
#include "ParticleBehavior.h"

class ParticleOrientation;
class SingularityInterrogator;
class ParticleShaderChain;
class ImplicitInterrogator;

class SingularityAdhesion : public ParticleBehavior {
	
private:
	
	// Particle orientations.
	ParticleOrientation* p_orient;
	// the target system's singularity chain
	ParticleShaderChain* sinChn;

	/// An implicit function to interrogate.
	ImplicitInterrogator* impInt;
	// a floater know its target particle system is the control
	Particles *target;
	
	/// parameters
	// population tolerance
	int popDifTolerance;
	// movement tolerance
	double movTolerance;
	// number of singularities tolerance
	int maxSingularPopulation;
	// number of time need to be stable tolerance
	int stableTolerance;
	// set to 1 to recaculate
	int recalculate;
	 // distance traveled more than this will be showed
	double distnaceTol;
	// angels in radient more than with will be showed
	double normalTol;

	/// flags to track when to start collect points
	// boolean to set if the system is stable
	bool isStable;
	// total movement
	double movement;

	/// distance traveled over last n iteration
	std::vector<double> distance;
	/// normal angle change in radient
	std::vector<double> normalChange;

	/// vectors to record the history
	/// singular positions
	std::vector<gmVector3> singularPoints;
	/// singular cumulation
	std::vector<int> singularIndex;
	/// normal angle change in radient
	std::vector<gmVector3> oldN;


public:

	MAKE_PARTICLESTUFF_NAME();

	void createSingularPS(Particles* ps);
	virtual void attachAttributes();
	
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

	int qlenpp();
	void getqpp(double *q, int i);
	void setqpp(double *q, int i);
	void qnamepp(char **qn);
	
	void particleAdded();
	
	void particleRemoved(unsigned int i);

	/// default constructor
	SingularityAdhesion(Particles *ps=NULL);
	
	/// destructor
	~SingularityAdhesion() {}

	/// Returns the particle orientation
	ParticleOrientation* getOrientation() { return p_orient; }

	/// record the old normal
	void applyForce();

	/// record the old normal
	void applyConstraint();
		
	/// record the new normal
	void cleanup();

	/// helper to reset enviroment to recalculate singularities
	void resetSingularities();

	// update the history size
	void updateHistorySize();

};

#endif
