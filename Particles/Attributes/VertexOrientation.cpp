/*
@file VertexOrientation.cpp
@author Jerry O. Talton III
@date 14 Jan. 2005
*/

#include "VertexOrientation.h"
#include "ParticleNormal.h"
#include "MeshInterrogator.h"
#include "Surface/OpenMesh/SurfaceMesh.h"

REGISTER_PARTICLESTUFF(VertexOrientation,"Attribute:VertexOrientation");

VertexOrientation::VertexOrientation(Particles *ps, const std::string& name)
	: ParticleNormal(ps,name)
{
	new Attached<MeshInterrogator>(this,&mi);
}

void VertexOrientation::attachAttributes()
{
	ParticleNormal::attachAttributes();

	if (!mi) return;

	SurfaceMesh *mesh = mi->getMesh();
	if (!mesh) return;

	unsigned int i = 0;
	for (SurfaceMesh::VertexIter v_it = mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it) {
		if (i >= ps->size()) {
			ps->addParticle();
		}
		SurfaceMesh::Normal nm(mesh->normal(v_it));
		
		gmVector3 n(nm[0], nm[1], nm[2]);
		
		setNormal(i, n);
		i++;
	}
}
