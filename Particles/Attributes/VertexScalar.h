/**	@file VertexScalar.h
	@auther John C. Hart
	@date 31 Dec. 2005 (I gotta get a life)
*/

#ifndef VERTEXSCALAR_H
#define VERTEXSCALAR_H

#include "Particles.h"
#include "ParticleScalar.h"
class MeshInterrogator;

/** VertexScalar is a ParticleScalar with added functionality the load and save
 * the scalars as a texture coordinate of the vertices of the SurfaceMesh
 * described by a MeshInterrogator.
 */
class VertexScalar : public ParticleScalar {
public:
	MAKE_PARTICLESTUFF_NAME();

	bool persist;
	bool store;
	MeshInterrogator *mi;

	VertexScalar(Particles *ps=NULL, const std::string& name=std::string("VertexScalar"));

	virtual void attachAttributes();
	virtual void prepare();

	/// Resets particle scalars to mesh vertex texcoord 1.
	virtual void loadScalars();

	/// Updates vertex texcoord 1 to current particle scalars.
	virtual void saveScalars();
};

class LoadVertexScalar : public PSParamButton::Callback
{
public:
	VertexScalar *vs;
	LoadVertexScalar(VertexScalar *me) {vs = me;}
	virtual void onbuttonpress() {vs->loadScalars();}
};

#endif
