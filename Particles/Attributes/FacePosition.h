/**	@file FacePosition.h
	@auther John C. Hart
	@date 4 Jan. 2005
*/

#ifndef FACEPOSITION_H
#define FACEPOSITION_H

#include "Particles.h"
#include "ParticlePosition.h"
class MeshInterrogator;

class FacePosition : public ParticlePosition {
	//added support for a random selection of particles on the faces
	//0 leads to the old behavior: one particle per face
	//the random selection is compatible with FaceOrientation
	int _maxNbParticles;
public:
	MAKE_PARTICLESTUFF_NAME();

	MeshInterrogator *mi;

	FacePosition(Particles *ps=NULL, const std::string &name=std::string("FacePosition"));

	virtual void attachAttributes();
};

#endif 


