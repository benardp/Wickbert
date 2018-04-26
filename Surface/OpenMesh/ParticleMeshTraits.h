/*
@file ParticleMeshTraits.h
@author wensu
@date 2003-09-06
This is for maintain particles on a mesh. Each vertex must be 
paired with a particle in the particle system.
*/

#ifndef PARTICLEMESHTRAITS_H
#define PARTICLEMESHTRAITS_H

// openmesh includes
#include "OpenMesh/Core/IO/BinaryHelper.hh"
#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"

#ifdef _MSC_VER
#  ifndef OM_STATIC_BUILD
#  define OM_STATIC_BUILD
#  endif
#  define INCLUDE_TEMPLATES
#  include <OpenMesh/Core/IO/IOInstances.hh>
#endif

struct ParticleMeshTraits : public OpenMesh::DefaultTraits
{
	typedef float Scaler;

	FaceAttributes(OpenMesh::Attributes::Normal|OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
	VertexAttributes(OpenMesh::Attributes::Normal|OpenMesh::Attributes::Status);

	FaceTraits
	{
		Scaler area;
	};
	EdgeTraits
	{
		Scaler length;
		//added by Matei Stroila to fix the Status problem:
		bool tag;
	};
};

typedef OpenMesh::TriMesh_ArrayKernelT<ParticleMeshTraits> ParticleOpenMesh;

#endif
