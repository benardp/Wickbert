#ifndef		_ASYMMETRIC_REPULSION
#define		_ASYMMETRIC_REPULSION

#include "ParticleBehavior.h"
#include "AdaptiveRepulsionData.h"

class ParticleLocality;
class ParticlePosition;
class ParticleVelocity;

class AsymmetricRepulsion : public ParticleBehavior
{
public:
	MAKE_PARTICLESTUFF_NAME();
	
	std::string source;			//	name of particles from which to apply the asymmetric repulsion
	//	the particles from which to apply the asymmetric repulsion
	Particles*	source_p;

	//	constructor
	AsymmetricRepulsion(Particles* ps = NULL);

	void attachAttributes();

	void applyForce();
//	void integrate();	//	use the parent class implementation
//	void cleanup();		//	use the parent class implementation

	//	factor for the asymmetric force

	AdaptiveRepulsionData* rep_data;

	ParticlePosition* position;
	ParticleVelocity* velocity;

	double factor;
};

#endif	//	_ASYMMETRIC_REPULSION
