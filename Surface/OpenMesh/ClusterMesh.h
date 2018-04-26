/*
@file ClusterMesh.h
@author wensu
@date 2003-11-23
This is for maintain clusters on a mesh.
*/

#ifndef CLUSTERMESH_H
#define CLUSTERMESH_H

#include "Surface/Surface.h"
#include "Surface/OpenMesh/ClusterMeshTraits.h"

class ClusterMesh: public Surface
{
public:

	ClusterOpenMesh mesh;

	/// bounding box corners
	ClusterOpenMesh::Point bbox_min , bbox_max;

	/// default constructor
	ClusterMesh();

	/// open a file and precomputer all attributes
	bool readFile(const char* filename);

	void assignRandomCluster(int n);

	// set cluster to -1 reassign the cost of each face
	void removeCluster();

	void computeBoundingBox();

	void computeFaceCenter();

	void computeDualEdgeLength();

	void computeEdgeAngle();

	void computeEdgeCost();

	/// naming convention
	static std::string registry_name; 
	virtual const std::string name()
	{
		return registry_name;
	}
	virtual void resetObjectName()
	{
		std::ostringstream objectNumberName; 
		objectNumberName << ":" << ++objectNumber; 
		objectName = name() + objectNumberName.str(); 
	}

};

#endif
