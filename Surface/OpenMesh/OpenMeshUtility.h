/*
@file OpenMeshUtility.h
@author wensu
@date 2003-10-30
This calculate edge length and face areas for openmesh.
*/

#ifndef OPENMESHUTILITY_H
#define OPENMESHUTILITY_H
#include <cmath>

template< typename MeshT >
class ComputeEdgeLength
{
public:
    ComputeEdgeLength(MeshT& m) : mesh(m) {}
    void operator()(typename MeshT::Edge& e)
    {
		

		
		typename MeshT::EdgeHandle eh=mesh.handle(e);
		
		typename MeshT::HalfedgeHandle he=mesh.halfedge_handle(eh,0);
		typename MeshT::Point to(mesh.point(mesh.to_vertex_handle(he)));
		typename MeshT::Point from(mesh.point(mesh.from_vertex_handle(he)));
		typename MeshT::Point edge(to-from);
		e.length=edge.length();
    }

private:
    MeshT& mesh;
};

// assumes it is a traingle mesh
template< typename MeshT >
class ComputeFaceArea
{
public:
    ComputeFaceArea(MeshT& m) : mesh(m) {}
    void operator()(typename MeshT::Face& f)
    {
		typename MeshT::FaceHandle fh=mesh.handle(f);

		typename MeshT::FaceVertexIter fv_it=mesh.fv_iter(fh);
		typename MeshT::VertexHandle va=fv_it.handle(); ++fv_it;
		typename MeshT::VertexHandle vb=fv_it.handle(); ++fv_it;
		typename MeshT::VertexHandle vc=fv_it.handle();
		typename MeshT::Point c(mesh.point(va)-mesh.point(vb));
		typename MeshT::Point b(mesh.point(va)-mesh.point(vc));
		// the absolute value
		f.area=(float)fabs((c%b).length()*0.5);
	}
private:
    MeshT& mesh;
};

template< typename MeshT >
class EdgeLengthGreaterThan
{
public:
    EdgeLengthGreaterThan(MeshT& m) : mesh(m) {}
	bool operator()(typename MeshT::EdgeHandle &e1, typename MeshT::EdgeHandle &e2)
	{
		return mesh.edge(e1).length > mesh.edge(e2).length;
	}

private:
    MeshT& mesh;
};

template< typename MeshT >
class HalfedgeLengthLessThan
{
public:
	HalfedgeLengthLessThan(MeshT &m) : mesh(m) {}

	bool operator()(typename MeshT::HalfedgeHandle &e1, typename MeshT::HalfedgeHandle &e2)
	{
		typename MeshT::EdgeHandle v1=mesh.edge_handle(e1);
		typename MeshT::EdgeHandle v2=mesh.edge_handle(e2);
		return mesh.edge(v1).length < mesh.edge(v2).length;
	}

private:
	MeshT& mesh;
};

template< typename MeshT >
class FaceAreaLessThan
{
public:
	FaceAreaLessThan(MeshT &m) : mesh(m) {}
	bool operator()(typename MeshT::FaceHandle &f1, typename MeshT::FaceHandle &f2)
	{
		return mesh.face(f1).area < mesh.face(f2).area;
	}
private:
	MeshT& mesh;
};

template< typename MeshT >
class FaceAreaGreaterThan
{
public:
	FaceAreaGreaterThan(MeshT &m) : mesh(m) {}
	bool operator()(typename MeshT::FaceHandle &f1, typename MeshT::FaceHandle &f2)
	{
		return mesh.face(f1).area > mesh.face(f2).area;
	}
private:
	MeshT& mesh;
};

#endif