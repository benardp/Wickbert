/*
@file ParticleCreation.cpp
@author Wen Su
*/

#include "ParticleMesh.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"

REGISTER_PARTICLESTUFF(ParticleMesh,"Attribute:ParticleMesh");

ParticleMesh::ParticleMesh(Particles *ps)
	: ParticlePosition(ps, std::string("ParticleMesh"), false)
{
	new PSParamString(this,&filename,"<empty>","filename","Mesh Source File",
		"Source file from which to input mesh.");
    new PSParamString(this, &outFileName,"", "outFileName", "Output Mesh File",
        "Output file for the generated mesh.");
    new PSParamButton(this, new WriteMeshCB(this), "Write Mesh", "Write Mesh Callback",
        "Writes the mesh to the file in Output Mesh File.");
}

void ParticleMesh::attachAttributes()
{
	if ((filename != "<empty>") && (filename != loaded_filename)) {
		loaded_filename = filename;
        ps->removeAll();
		mesh.clear();

		if (!OpenMesh::IO::read_mesh(mesh, filename.c_str()))
            return;

		for(ParticleOpenMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end();++v_it)
			ps->addParticle();
	}
}

/**
 * Overriding the particle position implementation. Nooping for now, it's difficult to
 * implement something as adding a particle has to consider the mesh data structure:
 * e.g. which vertices are neighbors to this new vertex?  Right now not enough
 * information is passed to particleAdded to communicate this.  Thus, we leave
 * everything up to the Behavior which decides something must be added.
 */
void ParticleMesh::particleAdded()
{
}

/**
 * Overriding the particle position implementation. Nooping for now, it's difficult to
 * implement something as removing a particle has to consider the mesh data structure:
 * e.g. how is the connectivity going to be updated?  Right now not enough
 * information is passed to particleRemoved to communicate this.  Thus, we leave
 * everything up to the Behavior which decides something must be removed.
 */
void ParticleMesh::particleRemoved(unsigned int i)
{
}

void ParticleMesh::createDiamond(double size)
{
	// generate vertices
	ParticleOpenMesh::VertexHandle vhandle[6];
	vhandle[0] = mesh.add_vertex(ParticleOpenMesh::Point( 0,  0,(ParticleMeshTraits::Scaler) -size));
	vhandle[1] = mesh.add_vertex(ParticleOpenMesh::Point( 0,  0, (ParticleMeshTraits::Scaler) size));
	vhandle[2] = mesh.add_vertex(ParticleOpenMesh::Point( 0,  (ParticleMeshTraits::Scaler) size,  0));
	vhandle[3] = mesh.add_vertex(ParticleOpenMesh::Point( 0, (ParticleMeshTraits::Scaler) -size,  0));
	vhandle[4] = mesh.add_vertex(ParticleOpenMesh::Point((ParticleMeshTraits::Scaler)  size,  0,  0));
	vhandle[5] = mesh.add_vertex(ParticleOpenMesh::Point((ParticleMeshTraits::Scaler) -size,  0,  0));

	// generate faces
	std::vector<ParticleOpenMesh::VertexHandle>  face_vhandles;
	face_vhandles.clear();
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[4]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[4]);
	face_vhandles.push_back(vhandle[3]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[3]);
	face_vhandles.push_back(vhandle[5]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[5]);
	face_vhandles.push_back(vhandle[2]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[5]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[5]);
	face_vhandles.push_back(vhandle[3]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[3]);
	face_vhandles.push_back(vhandle[4]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();
	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[4]);
	face_vhandles.push_back(vhandle[2]);
	mesh.add_face(face_vhandles);

	// update particles
	ps->addParticle();
	ps->addParticle();
	ps->addParticle();
	ps->addParticle();
	ps->addParticle();
	ps->addParticle();
}

void ParticleMesh::setPosition(unsigned int i, gmVector3 p)
{
	ParticleOpenMesh::VertexHandle vh(i);
	if (vh.is_valid())
		mesh.set_point(vh,ParticleOpenMesh::Point((ParticleMeshTraits::Scaler) p[0],(ParticleMeshTraits::Scaler) p[1],(ParticleMeshTraits::Scaler) p[2]));
}

gmVector3 ParticleMesh::getPosition(unsigned int i)
{
	ParticleOpenMesh::VertexHandle vh(i);
	if (vh.is_valid())
	{
		ParticleOpenMesh::Point p(mesh.point(vh));
		return gmVector3(p[0],p[1],p[2]);
	}
	return gmVector3();
}


void ParticleMesh::clear()
{
	mesh.clear();
}

void ParticleMesh::writeMesh()
{
    if (outFileName != "")
        OpenMesh::IO::write_mesh(mesh, outFileName.c_str());
}

ParticleMesh::WriteMeshCB::WriteMeshCB(ParticleMesh *mesh)
{
    _mesh = mesh;
}

void ParticleMesh::WriteMeshCB::onbuttonpress()
{
    _mesh->writeMesh();
}