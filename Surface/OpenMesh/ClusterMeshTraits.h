/*
@file ClusterMeshTraits.h
@author wensu
@date 2003-11-14
This is for maintain clusters on a mesh.
*/

#ifndef CLUSTERMESHTRAITS_H
#define CLUSTERMESHTRAITS_H

// openmesh includes
#include <OpenMesh/Core/IO/BinaryHelper.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#ifdef _MSC_VER
#  ifndef OM_STATIC_BUILD
#  define OM_STATIC_BUILD
#  endif
#  define INCLUDE_TEMPLATES
#  include <OpenMesh/Core/IO/IOInstances.hh>
#endif

struct ClusterMeshTraits : public OpenMesh::DefaultTraits
{

	FaceAttributes(OpenMesh::Attributes::Normal|OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
	VertexAttributes(OpenMesh::Attributes::Normal|OpenMesh::Attributes::Status);

	FaceTraits
	{
		Point center;
		float cost; // cost from cluster center to this face
		int cluster;
	};

	EdgeTraits
	{
		float dualEdgeLength;
		float angle;
		float cost; // cost of walking between two neighboring faces
	};
};

typedef OpenMesh::TriMesh_ArrayKernelT<ClusterMeshTraits> ClusterOpenMesh;

template< typename MeshT >
class FaceCostLessThan
{
public:
    FaceCostLessThan(MeshT& m) : mesh(m) {}
	bool operator()(typename MeshT::FaceHandle &f1, typename MeshT::FaceHandle &f2)
	{
		return mesh.face(f1).cost > mesh.face(f2).cost;
	}

private:
    MeshT& mesh;
};

#endif
