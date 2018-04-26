/**	@file VertexPosition.h
	@auther John C. Hart
	@date 4 Jan. 2005
*/

#ifndef OBJPOSITION_H
#define OBJPOSITION_H

#include "Particles.h"
#include "ParticlePosition.h"
class MeshInterrogator;

class OBJPosition : public ParticlePosition {
public:
	MAKE_PARTICLESTUFF_NAME();

	MeshInterrogator *mi;
	int numElements;
	int getNumElements(){ numElements = x.size(); return numElements;};
	std::string objFilename;

	OBJPosition(Particles *ps=NULL, const std::string& name=std::string("OBJPosition"));

	virtual void attachAttributes();
};

#endif

