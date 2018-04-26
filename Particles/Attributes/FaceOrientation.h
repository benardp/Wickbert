/**	@file FaceOrientation.h
	@auther Jerry O. Talton III
	@date 14 Jan. 2005
*/

#ifndef VERTEXORIENTATION_H
#define VERTEXORIENTATION_H

#include "Particles.h"
#include "ParticleNormal.h"
class MeshInterrogator;

class FaceOrientation : public ParticleNormal {
	//added support for a random selection of particles on the faces
	//0 leads to the old behavior: one particle per face
	//the random selection is compatible with FacePosition
	int _maxNbParticles;
public:
	MAKE_PARTICLESTUFF_NAME();

	MeshInterrogator *mi;

	FaceOrientation(Particles *ps=NULL, const std::string& name=std::string("FaceOrientation"));

	virtual void attachAttributes();
};

#endif

