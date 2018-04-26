/*
@file VertexPosition.cpp
@author John C. Hart
@date 4 Jan. 2005
*/

#include "VertexPosition.h"
#include "ParticlePosition.h"
#include "MeshInterrogator.h"
#include "Surface/OpenMesh/SurfaceMesh.h"

REGISTER_PARTICLESTUFF(VertexPosition,"Attribute:VertexPosition");

VertexPosition::VertexPosition(Particles *ps, const std::string& name)
	: ParticlePosition(ps,name)
{
	new PSParamBool(this,&persist,false,"persist","Persist",
		"If set, particle positions do not move from vertex locations.");
	new PSParamBool(this,&store,false,"store","Store",
		"If set, particle positions update vertex locations.");
	new PSParamButton(this, new LoadVertexPosition(this),"load","Reset Positions",
		"Reset particle positions to mesh vertex locations.");

	// new Attached<ParticlePosition>(this,&position);
	new Attached<MeshInterrogator>(this,&mi);
}

void VertexPosition::attachAttributes()
{
	ParticlePosition::attachAttributes();
	loadPositions();
}

void VertexPosition::prepare()
{
	if (persist) loadPositions();
	if (store) savePositions();
}

void VertexPosition::loadPositions()
{
	if (!mi) return;

	SurfaceMesh *mesh = mi->getMesh();
	if (!mesh) return;

	unsigned int i = 0;
	for (SurfaceMesh::VertexIter v_it = mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it) {
		if (i >= ps->size()) {
			ps->addParticle();
		}
		SurfaceMesh::Point p(mesh->point(v_it));
		x[i][0] = p[0];
		x[i][1] = p[1];
		x[i][2] = p[2];

		i++;
	}
	for (unsigned int j = ps->size() - 1; j >= i; j--)
		ps->removeParticle(j);
}

void VertexPosition::savePositions()
{
	if (!mi) return;

	SurfaceMesh *mesh = mi->getMesh();
	if (!mesh) return;

	SurfaceMesh::VertexIter v_it = mesh->vertices_begin(); 
	for (unsigned int i = 0; i < ps->size(); i++) {
		if (v_it == mesh->vertices_end()) {
			mesh->add_vertex(MyMesh::Point(	(SurfaceMesh::Scalar)x[i][0]
											,(SurfaceMesh::Scalar)x[i][1]
											,(SurfaceMesh::Scalar)x[i][2]));
		} else {
			mesh->set_point(v_it, MyMesh::Point( (SurfaceMesh::Scalar)x[i][0]
												,(SurfaceMesh::Scalar)x[i][1]
												,(SurfaceMesh::Scalar)x[i][2]));
			++v_it;
		}
	}
	while (v_it != mesh->vertices_end()) {
		SurfaceMesh::VertexIter extra = v_it;
		++v_it;	// could do this above but may be confusing
		mesh->delete_vertex(extra);
	}
}
