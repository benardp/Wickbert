/**
 * @file MeshShape.cpp
 * @date 2003-09-07
 * @author Wen Su
 */

#include "MeshShape.h"
#include "ParticleMesh.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "Surface/OpenMesh/OpenMeshUtility.h"
#include <algorithm>

REGISTER_PARTICLESTUFF(MeshShape,"Behavior:MeshShape");

MeshShape::MeshShape(Particles *ps)
:ParticleBehavior(ps, std::string("MeshShape"))
{
#if 0
	edgeLengthScaler=0.5;
	targetTriangleArea=0.2;
	longestWeightTolerance=1.3;
	shortestWeightTolerance=0.7;
	splitStepLimit=0.05;
	collapseStepLimit=0.05;
#else
	new PSParamDouble(this,&edgeLengthScaler,0.5,
		"edgescale","edge length scalar","Velocity created by differences in edge lengths");
	new PSParamDouble(this,&targetTriangleArea,0.2,
		"triarea","target triangle area","Goal triangle size");
	new PSParamDouble(this,&longestWeightTolerance,1.3,
		"longwghttol","longest weight tolerance","High end of edge weights");
	new PSParamDouble(this,&shortestWeightTolerance,1.3,
		"shortwghttol","shortest weight tolerance","Low end of edge weights");
	new PSParamDouble(this,&splitStepLimit,0.05,
		"splitsteplimit","split step limit","Portion of worst edges on which to focus");
	new PSParamDouble(this,&collapseStepLimit,0.05,
		"collapsesteplimit","collapse step limit","Portion of all edges to collapse");
#endif
}

/// each vertex find its longest edge move to it if it differs greatly
/// each vertex find its shortest edge move to it if it differs greatly
void MeshShape::applyForce()
{
	// position must be a ParticleMesh.
	ParticleMesh *pm=dynamic_cast<ParticleMesh *>(position);
	if (!pm)
		return;
	ParticleOpenMesh &mesh=pm->mesh;

	for(ParticleOpenMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end();++v_it)
	{
		// find average length
		double totalLength=0;
		// find longest edge
		double longest=0;
		double shortest=1e50;
		double edgeLength;
		double weight;
		int valence=0;
		ParticleOpenMesh::Point longForceVector;
		ParticleOpenMesh::Point shortForceVector;
		ParticleOpenMesh::Point centerVertex(mesh.point(v_it.handle()));
		for (ParticleOpenMesh::VertexVertexIter vv_it=mesh.vv_iter(v_it.handle()); vv_it;++vv_it)
		{
			ParticleOpenMesh::Point edgeVector(mesh.point(vv_it.handle())-centerVertex);
			edgeLength=edgeVector.length();
			totalLength+=edgeLength;
			if (edgeLength>longest)
			{
				longest=edgeLength;
				longForceVector=edgeVector;
			}
			if (edgeLength<shortest)
			{
				shortest=edgeLength;
				shortForceVector=edgeVector;
			}
			++valence;
		}
		//std::cout << "check valence" << valence << std::endl;

		weight=valence*longest/totalLength;
		// move to that direction
		if (weight>longestWeightTolerance)
		{
			velocity->v[v_it.handle().idx()]+=edgeLengthScaler*
				gmVector3(longForceVector[0],longForceVector[1],longForceVector[2]);
		}

		weight=valence*shortest/totalLength;
		// move to that direction
		if (weight<shortestWeightTolerance)
		{
			velocity->v[v_it.handle().idx()]-=edgeLengthScaler*
				gmVector3(shortForceVector[0],shortForceVector[1],shortForceVector[2]);
		}
	}
}

/// triangle edge flips
void MeshShape::applyConstraint()
{
	// position must be a ParticleMesh.
	ParticleMesh *pm=dynamic_cast<ParticleMesh *>(position);
	if (!pm)
		return;
	ParticleOpenMesh &mesh=pm->mesh;

	// computer the edge lengths 
	std::for_each(mesh.edges_begin(), mesh.edges_end(), ComputeEdgeLength<ParticleOpenMesh>(mesh));

	std::vector<ParticleOpenMesh::HalfedgeHandle> flipEdgeList;
	for (ParticleOpenMesh::EdgeIter e_it=mesh.edges_begin(); e_it!=mesh.edges_end();++e_it)
	{
		// if the dual edge is shorter than this edge do a flip.
		ParticleOpenMesh::HalfedgeHandle he=mesh.halfedge_handle(e_it.handle(),0);
		ParticleOpenMesh::HalfedgeHandle ohe=mesh.halfedge_handle(e_it.handle(),1);
		ParticleOpenMesh::VertexHandle top=mesh.to_vertex_handle(mesh.next_halfedge_handle(he));
		ParticleOpenMesh::VertexHandle down=mesh.to_vertex_handle(mesh.next_halfedge_handle(ohe));

		ParticleOpenMesh::Point dualEdgeLength=mesh.point(top)-mesh.point(down);
		float edgeLength=mesh.edge(e_it.handle()).length;
		float diff=edgeLength - dualEdgeLength.length();
		// tolerance prevent constant flipping
		if (diff/edgeLength > 0.05)
		{
			//std::cout << "edge length diff " << diff << std::endl;
			//std::cout << "edge length " << edgeLength.length() << std::endl;
			flipEdgeList.push_back(he);
		}
	}

	//if (flipEdgeList.size()>0)
	//	std::cout << "flipEdgeList size: " <<flipEdgeList.size() << std::endl << "list: " ;
	for(unsigned int i=0;i<flipEdgeList.size();++i)
	{
		//std::cout << flipEdgeList[i] << ", ";
		// this uses OpenMesh function
		if (pm->mesh.is_flip_ok(pm->mesh.edge_handle(flipEdgeList[i])))
			pm->mesh.flip(pm->mesh.edge_handle(flipEdgeList[i]));
	}
	//if (flipEdgeList.size()>0)
	//	std::cout << std::endl;
}


///Depending target area, it splits or collapses
void MeshShape::cleanup()
{
	// position must be a ParticleMesh.
	ParticleMesh *pm=dynamic_cast<ParticleMesh *>(position);
	if (!pm)
		return;
	ParticleOpenMesh &mesh=pm->mesh;

	// computer the edge lengths 
	// sorting depend on this values
	std::for_each(mesh.faces_begin(), mesh.faces_end(), ComputeFaceArea<ParticleOpenMesh>(mesh));

	// split and collapse limit base on target area
	double splitLimit=1.5*targetTriangleArea;
	double collapseLimit=0.5*targetTriangleArea;

	// subdived if larger than target area.
	// cannot collapse two neighboring faces
	std::vector<ParticleOpenMesh::FaceHandle> split;
	std::vector<ParticleOpenMesh::FaceHandle> collapse;
	for (ParticleOpenMesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end(); ++f_it)
	{
		if (f_it->area>splitLimit)
			split.push_back(f_it.handle());
		else if (f_it->area<collapseLimit)
			collapse.push_back(f_it.handle());
	}

	// sort only top splitStepLimit of the number
	unsigned int splitSortSize=splitStepLimit*split.size()+1;
	// 0.05*10 = 0.5, so still should split at least 1.
	splitSortSize=(split.size()<splitSortSize?split.size():splitSortSize);
	std::partial_sort(split.begin(), split.begin()+splitSortSize, split.end(), FaceAreaGreaterThan<ParticleOpenMesh>(mesh));
	// only split, hopefully edge flips will take care of the bad shapes.
	for(unsigned int i=0;i<splitSortSize;i++)
	{
		//std::cout << "split face : "<< i <<" area: " << mesh.face(split[i]).area <<std::endl;
		ParticleOpenMesh::FaceVertexIter fv_it=mesh.fv_iter(split[i]);
		ParticleOpenMesh::VertexHandle va=fv_it.handle(); ++fv_it;
		ParticleOpenMesh::VertexHandle vb=fv_it.handle(); ++fv_it;
		ParticleOpenMesh::VertexHandle vc=fv_it.handle();
		ParticleOpenMesh::Point center(mesh.point(va)+mesh.point(vb)+mesh.point(vc));
		center=center*0.333333333333333333333f;
		mesh.vertex_split(center,va,vb,vc);
	 	ps->addParticle();
	}

	// sort only top splitStepLimit of the number
	unsigned int collapseSortSize=collapseStepLimit*collapse.size()+1;
	collapseSortSize=(collapse.size()<collapseSortSize?collapse.size():collapseSortSize);
	std::partial_sort(collapse.begin(), collapse.begin()+collapseSortSize, collapse.end(), FaceAreaLessThan<ParticleOpenMesh>(mesh));

	// only split, hopefully edge flips will take care of the bad shapes.
	for(unsigned int i=0;i<collapseSortSize;i++)
	{
		if (collapse[i].idx()<0 || collapse[i].idx()>mesh.n_faces())
			continue;
		ParticleOpenMesh::HalfedgeHandle he=mesh.halfedge_handle(collapse[i]);
		if (he.idx()<0 || he.idx()>mesh.n_edges())
			continue;
		// this uses OpenMesh function
		if (pm->mesh.is_collapse_ok(he))
		{
			//std::cout << "collapse face : "<< i <<" area: " << mesh.face(collapse[i]).area <<std::endl;
			ParticleOpenMesh::VertexHandle nvh=pm->mesh.from_vertex_handle(he);
			pm->mesh.collapse(he);
			pm->mesh.garbage_collection();
			// Halfedge collapse: collapse the from-vertex of halfedge _heh into its to-vertex. 
	 		ps->removeParticle(nvh.idx());
		}
	}

}

int MeshShape::qlen()
{
	return 6;
}
void MeshShape::getq(double *q)
{
	q[0] = edgeLengthScaler;
	q[1] = targetTriangleArea;
	q[2] = longestWeightTolerance;
	q[3] = shortestWeightTolerance;
	q[4] = splitStepLimit;
	q[5] = collapseStepLimit;
}
void MeshShape::setq(double *q)
{
	edgeLengthScaler = q[0];
	targetTriangleArea = q[1];
	longestWeightTolerance = q[2];
	shortestWeightTolerance = q[3];
	splitStepLimit = q[4];
	collapseStepLimit =q[5];
}
void MeshShape::qname(char **qn)
{
	qn[0] = "edgeLengthScaler";
	qn[1] = "targetTriangleArea";
	qn[2] = "longestWeightTolerance";
	qn[3] = "shortestWeightTolerance";
	qn[4] = "splitStepLimit (percentage)";
	qn[5] = "collapseStepLimit (percentage)";
}
