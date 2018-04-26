/**	@file VertexOrientation.h
	@auther Jerry O. Talton III
	@date 14 Jan. 2005
*/

#ifndef VERTEXORIENTATION_H
#define VERTEXORIENTATION_H

#include "Particles.h"
#include "ParticleNormal.h"
class MeshInterrogator;

class VertexOrientation : public ParticleNormal {
public:
	MAKE_PARTICLESTUFF_NAME();

	MeshInterrogator *mi;

	VertexOrientation(Particles *ps=NULL, const std::string& name=std::string("VertexOrientation"));

	virtual void attachAttributes();
};

#endif

