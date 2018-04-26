/*
@file ClusterMesh.cpp
@author wensu
@date 2003-11-23
This is for maintain clusters on a mesh.
*/

#include "ClusterMesh.h"
#include <cmath>
#include <limits>

std::string ClusterMesh::registry_name = "ClusterMesh";

ClusterMesh::ClusterMesh()
:Surface()
{
}

bool ClusterMesh::readFile(const char* filename)
{
	if (filename==NULL)
		return false;
	mesh.clear();
	bool result=OpenMesh::IO::read_mesh(mesh, filename);

	// some debug information
	std::cout << "faces: "<< mesh.n_faces()<<std::endl;
	std::cout << "vertices: "<< mesh.n_vertices()<<std::endl;
	std::cout << "edges: "<< mesh.n_edges()<<std::endl;

	// precomputing step
	removeCluster();
	mesh.update_face_normals();
	mesh.update_vertex_normals();
	computeBoundingBox();
	computeFaceCenter();
	computeDualEdgeLength();
	computeEdgeAngle();
	computeEdgeCost();
	return result;
}

void ClusterMesh::assignRandomCluster(int n)
{
	int base=mesh.n_faces()/n;
	int i=0;
	for (ClusterOpenMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
	{
		mesh.data(*fit).cluster=i++/base;
	}
}

void ClusterMesh::removeCluster()
{
	for (ClusterOpenMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
	{
		mesh.data(*fit).cost=std::numeric_limits<float>::max();
		mesh.data(*fit).cluster=-1;
	}
}

void ClusterMesh::computeBoundingBox()
{
	ClusterOpenMesh::VertexIter vit = mesh.vertices_begin(), vend = mesh.vertices_end();

	if (vit == vend)
		return;

	
	ClusterOpenMesh::Point firstPoint = mesh.point(vit);
	bbox_min = bbox_max = firstPoint;

	for ( ; vit != vend; ++vit)
	{		
		
		ClusterOpenMesh::Point p = mesh.point(vit);
		for (int i=0; i < 3; i++)
		{
			bbox_min[i] = bbox_min[i]<p[i]?bbox_min[i]:p[i];
			bbox_max[i] = bbox_max[i]>p[i]?bbox_min[i]:p[i];
		}
	}
}

void ClusterMesh::computeFaceCenter()
{
	for (ClusterOpenMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
	{
		ClusterOpenMesh::Point center(0,0,0);
		for (ClusterOpenMesh::FaceVertexIter fv_it=mesh.fv_iter(fit.handle()); fv_it ; ++fv_it)
		{
			center+=mesh.point(fv_it);
		}
		mesh.data(*fit).center=center/3;
	}
}
	
void ClusterMesh::computeDualEdgeLength()
{
	// boundaries edge is 0
	for (ClusterOpenMesh::EdgeIter eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit)
	{
		if (mesh.is_boundary(eit.handle()))
		{
			mesh.data(*eit).dualEdgeLength=0;
			continue;
		}
		ClusterOpenMesh::HalfedgeHandle he=mesh.halfedge_handle(eit.handle(), 0);
		ClusterOpenMesh::HalfedgeHandle ohe=mesh.halfedge_handle(eit.handle(), 1);
		ClusterOpenMesh::Point from(mesh.data(mesh.face_handle(he)).center);
		ClusterOpenMesh::Point to(mesh.data(mesh.face_handle(ohe)).center);
		ClusterOpenMesh::Point distance(from-to);
		mesh.data(*eit).dualEdgeLength=distance.length();
	}
}

void ClusterMesh::computeEdgeAngle()
{
	for (ClusterOpenMesh::EdgeIter eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit)
	{
		if (mesh.is_boundary(eit.handle()))
		{
			mesh.data(*eit).dualEdgeLength=0;
			continue;
		}
		ClusterOpenMesh::HalfedgeHandle he=mesh.halfedge_handle(eit.handle(), 0);
		ClusterOpenMesh::HalfedgeHandle ohe=mesh.halfedge_handle(eit.handle(), 1);
		
		ClusterOpenMesh::Point from(mesh.normal(mesh.face_handle(he)));
	
		ClusterOpenMesh::Point to(mesh.normal(mesh.face_handle(ohe)));
		// all normal are unit vectors so n.v=cos alpha
		mesh.data(*eit).angle=(from[0]*to[0]+from[1]*to[1]+from[2]*to[2]);
	}
}

void ClusterMesh::computeEdgeCost()
{
	for (ClusterOpenMesh::EdgeIter eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit)
	{
		mesh.data(*eit).cost=(1-mesh.data(*eit).angle)*(mesh.data(*eit).dualEdgeLength);
		//std::cout << eit->angle << " " << eit->dualEdgeLength << " " << eit->cost <<std::endl;
	}
}
