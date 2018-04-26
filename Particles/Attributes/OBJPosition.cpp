/*
@file OBJPosition.cpp
@author John C. Hart
@date 4 Jan. 2005
*/

#include "OBJPosition.h"
#include "ParticlePosition.h"
#include "MeshInterrogator.h"
#include "Surface/OpenMesh/SurfaceMesh.h"

REGISTER_PARTICLESTUFF(OBJPosition,"Attribute:OBJPosition");

OBJPosition::OBJPosition(Particles *ps, const std::string& name)
	: ParticlePosition(ps,name)
{
	new PSParamInt(this,&numElements,0,"nmElements","Number of elements", "Number of elements in this attribute");
	new PSParamString(this,&objFilename,"rbfDefaultCenters.obj","objFileNm","obj values filename","the name of the file that holds the vertex points of obj files");
	
	// new Attached<ParticlePosition>(this,&position);
	new Attached<MeshInterrogator>(this,&mi);
}

void OBJPosition::attachAttributes()
{
	ParticlePosition::attachAttributes();

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
}