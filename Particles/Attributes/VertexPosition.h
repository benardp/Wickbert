/**	@file VertexPosition.h
	@auther John C. Hart
	@date 4 Jan. 2005
*/

#ifndef VERTEXPOSITION_H
#define VERTEXPOSITION_H

#include "Particles.h"
#include "ParticlePosition.h"
class MeshInterrogator;

/** VertexPosition is a ParticlePosition with added functionality the load and save
 * positions as the vertices of the SurfaceMesh described by a MeshInterrogator.
 */
class VertexPosition : public ParticlePosition {
public:
	MAKE_PARTICLESTUFF_NAME();

	bool persist;
	bool store;
	MeshInterrogator *mi;

	VertexPosition(Particles *ps=NULL, const std::string& name=std::string("VertexPosition"));

	virtual void attachAttributes();
	virtual void prepare();

	/// Resets particle positions to mesh vertex locations.
	virtual void loadPositions();

	/// Updates vertex locations to current particle positions.
	virtual void savePositions();
};

class LoadVertexPosition : public PSParamButton::Callback
{
public:
	VertexPosition *vp;
	LoadVertexPosition(VertexPosition *me) {vp = me;}
	virtual void onbuttonpress() {vp->loadPositions();}
};

#endif


