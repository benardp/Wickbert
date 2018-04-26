/*
@file VertexScalar.cpp
@author John C. Hart
@date 31 Dec. 2005 (I gotta get a life)
*/

#include "VertexScalar.h"
#include "MeshInterrogator.h"
#include "Surface/OpenMesh/SurfaceMesh.h"

REGISTER_PARTICLESTUFF(VertexScalar,"Attribute:VertexScalar");

VertexScalar::VertexScalar(Particles *ps, const std::string& name)
	: ParticleScalar(ps,name)
{
	new PSParamBool(this,&persist,false,"persist","Persist",
		"If set, particle scalars do not move from vertex texcoord.");
	new PSParamBool(this,&store,false,"store","Store",
		"If set, particle scalars update vertex texcoord.");
	new PSParamButton(this, new LoadVertexScalar(this),"load","Reset Scalars",
		"Reset particle scalars to mesh vertex texcoord.");

	// new Attached<ParticlePosition>(this,&position);
	new Attached<MeshInterrogator>(this,&mi);
}

void VertexScalar::attachAttributes()
{
	ParticleScalar::attachAttributes();
	loadScalars();
}

void VertexScalar::prepare()
{
	if (persist) loadScalars();
	if (store) saveScalars();
}

void VertexScalar::loadScalars()
{
	if (!mi) return;

	SurfaceMesh *mesh = mi->getMesh();
	if (!mesh) return;

	unsigned int i = 0;
	for (SurfaceMesh::VertexIter v_it = mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it) {
		if (i >= ps->size()) {
			ps->addParticle();
		}
		SurfaceMesh::TexCoord2D tc = mesh->texcoord2D(v_it);
		
		setScalar(i,tc[0]);

		i++;
	}
	for (unsigned int j = ps->size() - 1; j >= i; j--)
		ps->removeParticle(j);
}

void VertexScalar::saveScalars()
{
	if (!mi) return;

	SurfaceMesh *mesh = mi->getMesh();
	if (!mesh) return;

	SurfaceMesh::VertexIter v_it = mesh->vertices_begin(); 
	for (unsigned int i = 0; i < ps->size(); i++) {
		if (v_it == mesh->vertices_end()) {
			mesh->add_vertex(MyMesh::Point());
		} else {
			SurfaceMesh::TexCoord2D tc = mesh->texcoord2D(v_it);
			tc[0] = (SurfaceMesh::Scalar) getScalar(i);
			mesh->set_texcoord2D(v_it, tc);
			++v_it;
		}
	}
	while (v_it != mesh->vertices_end()) {
		SurfaceMesh::VertexIter extra = v_it;
		++v_it;	// could do this above but may be confusing
		mesh->delete_vertex(extra);
	}
}
