//@file SPHConstraint.h
//@author Yuan Zhou

#ifndef SPHCONSTRAINT_H
#define SPHCONSTRAINT_H

#include "Particles.h"
#include "ParticleBehavior.h"
#include "ParticleBoundingBox.h"

class SPHConstraint: public ParticleBehavior
{
private:
	ParticleBoundingBox *bbox;

public:
	MAKE_PARTICLESTUFF_NAME();

	char *tip();

	/// Create a bounding box.
	SPHConstraint(Particles *ps=NULL, const std::string& name=std::string("SPHConstraint"));

	/// attach ParticleBoundingBox
	void attachAttributes();

	/// Constrains particles so that they do not leave the box.
	void applyConstraint();
};

#endif
