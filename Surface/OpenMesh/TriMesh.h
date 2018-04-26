/*
@file TriMesh.h
@author wensu
@date 2003-05-01
This allows a mesh to be used as an implicit surface. Using the offset of this is useful.
To use openmesh, please place it as the same dir as surface.
Use the default name OpenMesh_0.11.1 as the dir name.
Then surface should find both the include and library files OK.
*/


#ifndef TRIMESH_H
#define TRIMESH_H

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

#include <vector>
#include <limits>

struct SurfaceTraits : public OpenMesh::DefaultTraits
{
	FaceAttributes(OpenMesh::Attributes::Normal);
	//HalfedgeAttributes(OpenMesh::Attributes::Status);
	//EdgeAttributes(OpenMesh::Attributes::Status);
	VertexAttributes(OpenMesh::Attributes::Normal);

	FaceTraits
	{
		Point center;
	};

};

typedef OpenMesh::TriMesh_ArrayKernelT<SurfaceTraits> MyMesh;

#include "Surface/Implicit/Implicit.h"

class TriMesh: public Implicit
{
public:

	// default constructor
	TriMesh();

	MyMesh mesh;

	// bounding box
	MyMesh::Point bbox_min, bbox_max;

	// Read a file
	bool readFile(const char* filename=NULL);

	// Write a file
	bool writeFile(const char* filename=NULL);
	
	void computeFaceCenter();
	
	void computeBoundingBox();

	// return the distance from point p to line segment ab
	float distancePointToLine(const MyMesh::Point &p, const MyMesh::Point &a, const MyMesh::Point &b);

	// project Point p to plane formed by a b c
	MyMesh::Point projectPointToPlane(const MyMesh::Point &p, const MyMesh::Point &a, const MyMesh::Point &b, const MyMesh::Point &c);

	// project Point p to plane formed by point a and normal n
	MyMesh::Point projectPointToPlane(const MyMesh::Point &p, const MyMesh::Point &a, const MyMesh::Point &n);

	// not normalized normal of plane formed by a,b,c
	MyMesh::Point normalOf(const MyMesh::Point &a, const MyMesh::Point &b, const MyMesh::Point &c);

	bool isInTriangle(const MyMesh::Point &p, MyMesh::FaceHandle f);

	bool isInTriangle(MyMesh::VertexHandle p, MyMesh::VertexHandle a, MyMesh::VertexHandle b, MyMesh::VertexHandle c);

	bool isInTriangle(const MyMesh::Point &p, const MyMesh::Point &a, const MyMesh::Point &b, const MyMesh::Point &c);

	bool isOnSameSide(const MyMesh::Point &p1, const MyMesh::Point &p2, const MyMesh::Point &v1, const MyMesh::Point &v2);

	// closest distance to a point
	double distancePointToMesh(const MyMesh::Point &p);

	// distance to triangle
	double distancePointToTriangle(const MyMesh::Point &p, MyMesh::FaceHandle fh);

	MyMesh::VertexHandle closestVertex(const MyMesh::Point &p);

	MyMesh::FaceHandle closestFace(const MyMesh::Point &p);

    virtual double proc(const gmVector3 & x);

	virtual Intervald proc(const Box<double>& x);
	
	virtual gmVector3 normal(const gmVector3 & x);

	double exactDistanceToMesh(const MyMesh::Point &p);

//	virtual gmVector3 grad(const gmVector3 & x);

};

#endif
