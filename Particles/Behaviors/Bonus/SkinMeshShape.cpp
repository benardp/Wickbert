/**
 * @file SkinSkinMeshShape.cpp
 * @date 2003-09-07
 * @author Wen Su
 */

#include "Surface/OpenMesh/OpenMeshUtility.h"
#include "SkinMeshShape.h"
#include "ParticleVelocity.h"
#include "ImplicitInterrogator.h"
#include <cmath>
#include <algorithm>

REGISTER_PARTICLESTUFF(SkinMeshShape,"Behavior:SkinMeshShape");

SkinMeshShape::SkinMeshShape(Particles *ps)
:ParticleBehavior(ps, std::string("SkinMeshShape"))
{
	new PSParamDouble(this,&targetLength,0.8,
		"targetLength","target edge length scalar","Desired edge length for elements of the mesh.");
	new PSParamInt(this,&allowSplit, 5,
		"allowSplit","splits per iteration","Max number of edge splits per iteration.");
	new PSParamInt(this,&allowCollapse, 5,
		"allowCollapse","collapses per iteration","Max number of edge collapses per iteration.");
	new PSParamBool(this,&allowFlip,true,
		"allowFlip","allow edge flips","Desired edge length for elements of the mesh.");

    new Attached<ImplicitInterrogator>(this,&impInt);
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
}

void SkinMeshShape::attachAttributes()
{
	ParticleBehavior::attachAttributes();
}

/// new position=a*current position+b*neighbor+c*push to surface
void SkinMeshShape::applyForce()
{
	// position must be a ParticleMesh.
	ParticleMesh *pm=dynamic_cast<ParticleMesh *>(position);
	if (!pm)
		return;
	ParticleOpenMesh &mesh=pm->mesh;

	// for each particle
	for(ParticleOpenMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end();++v_it)
	{
		unsigned int valence=0;
		ParticleOpenMesh::VertexHandle vh(v_it.handle());
		ParticleOpenMesh::Point currentPosition=mesh.point(vh);
		ParticleOpenMesh::Point neighborPosition(0,0,0);
		for(ParticleOpenMesh::VertexVertexIter vv_it=mesh.vv_iter(vh); vv_it; ++vv_it)
		{
			++valence;
			// average position
			neighborPosition+=mesh.point(vv_it);
		}

		double m=0.0;	// estimate the smoothness
		mesh.update_face_normals();
		for (ParticleOpenMesh::VertexEdgeIter ve_it=mesh.ve_iter(vh); ve_it; ++ve_it)
		{
			// edge e1 to e2 has two faces
			ParticleOpenMesh::HalfedgeHandle e1(mesh.halfedge_handle(ve_it.handle(),0));
			ParticleOpenMesh::HalfedgeHandle e2(mesh.halfedge_handle(ve_it.handle(),1));
			ParticleOpenMesh::FaceHandle f1(mesh.face_handle(e1));
			ParticleOpenMesh::FaceHandle f2(mesh.face_handle(e2));
			
			ParticleOpenMesh::Point n1(mesh.normal(f1));
			ParticleOpenMesh::Point n2(mesh.normal(f2));
			
			m+=n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
		}
		neighborPosition=neighborPosition/((float)valence);
		m=m/valence;
		double a=0.3;
		double sq=(0.5*(m+1));
		double bScale=0.3*0.7*std::sqrt(1-sq*sq);
		double b=(1-a)*bScale;
		double c=1-a-b;
		double d=impInt->getImplicit()->proc(gmVector3(currentPosition[0],currentPosition[1],currentPosition[2]));
		gmVector3 normal(impInt->getImplicit()->normal(gmVector3(currentPosition[0],currentPosition[1],currentPosition[2])));
		ParticleOpenMesh::Point dNormal(	 (ParticleOpenMesh::Scalar)(normal[0]*d)
											,(ParticleOpenMesh::Scalar)(normal[1]*d)
											,(ParticleOpenMesh::Scalar)(normal[2]*d)
										);
		ParticleOpenMesh::Point surfacePosition(currentPosition-dNormal);
		ParticleOpenMesh::Point newPosition(	currentPosition*((ParticleOpenMesh::Scalar)a)
												+neighborPosition*((ParticleOpenMesh::Scalar)b)
												+surfacePosition*((ParticleOpenMesh::Scalar)c));
		ParticleOpenMesh::Point newVelocity(newPosition-currentPosition);
		//mesh.set_point(vh,newPosition);
		velocity->v[vh.idx()]+=gmVector3(newVelocity[0],newVelocity[1],newVelocity[2]);
	}
}

void SkinMeshShape::edgeCollapseAndSplit(ParticleOpenMesh &mesh)
{
	// clear status bits
	//for(ParticleOpenMesh::EdgeIter e_it=mesh.edges_begin(); e_it!=mesh.edges_end();++e_it)
	//	e_it->set_tagged(false);
	for(ParticleOpenMesh::EdgeIter e_it=mesh.edges_begin(); e_it!=mesh.edges_end();++e_it)
		e_it->tag = false;

	// computer the lengths and faces
	// sorting depend on this values
	std::for_each(mesh.edges_begin(), mesh.edges_end(), ComputeEdgeLength<ParticleOpenMesh>(mesh));

	double maxLength=1.5*targetLength;
	double minLength=0.5*targetLength;

	std::vector<ParticleOpenMesh::EdgeHandle> split;
	std::vector<ParticleOpenMesh::HalfedgeHandle> collapse;

	for(ParticleOpenMesh::EdgeIter e_it=mesh.edges_begin(); e_it!=mesh.edges_end();++e_it)
	{
		// do not flip on two neighboring edges
		if (e_it->tag)
			continue;

		ParticleOpenMesh::HalfedgeHandle he=mesh.halfedge_handle(e_it.handle(),0);
		ParticleOpenMesh::VertexHandle vha=mesh.from_vertex_handle(he);
		ParticleOpenMesh::VertexHandle vhb=mesh.to_vertex_handle(he);
		ParticleOpenMesh::Point edge(mesh.point(vha)-mesh.point(vhb));
		double length=edge.length();
		// split if greater than 1.5 of target length
		if (length > maxLength)
		{
			ParticleOpenMesh::Point midPoint(mesh.point(vha));
			midPoint+=mesh.point(vhb);
			midPoint=midPoint*0.5;
			split.push_back(e_it.handle());
		}
		// collapse if less than 0.5 of target length
		else if (length < minLength)
		{
			collapse.push_back(he);
		}
	}

	// only does top allowSplit
	unsigned int splitSortSize=((unsigned int)allowSplit<split.size()?(unsigned int)allowSplit:split.size());
	// split the longest n edges
	std::partial_sort(split.begin(), split.begin()+splitSortSize, split.end(), EdgeLengthGreaterThan<ParticleOpenMesh>(mesh));
	// check length
	//std::cout << "sort by length :" << split.size() << std::endl;
	for(unsigned int i=0;i<splitSortSize;++i)
	{
		ParticleOpenMesh::VertexHandle v1=mesh.from_vertex_handle(mesh.halfedge_handle(split[i],0));
		ParticleOpenMesh::VertexHandle v2=mesh.from_vertex_handle(mesh.halfedge_handle(split[i],1));
		ParticleOpenMesh::Point edge(mesh.point(v1)-mesh.point(v2));
		ParticleOpenMesh::Point center(mesh.point(v1)+mesh.point(v2));
		center=center*0.5;
		//std::cout << edge.length() << std::endl;
		mesh.split(split[i],mesh.add_vertex(center));
		ps->addParticle();
	}

	unsigned int collapseSortSize=((unsigned int)allowCollapse<collapse.size()?(unsigned int)allowCollapse:collapse.size());
	// collapse the shortest n edges
	std::partial_sort(collapse.begin(), collapse.begin()+collapseSortSize, collapse.end(), HalfedgeLengthLessThan<ParticleOpenMesh>(mesh));
	for(unsigned int i=0;i<collapseSortSize;++i)
	{
		if (mesh.is_collapse_ok(collapse[i]))
		{
			ParticleOpenMesh::VertexHandle nvh=mesh.from_vertex_handle(collapse[i]);
			mesh.collapse(collapse[i]);
			mesh.garbage_collection();
			// Halfedge collapse: collapse the from-vertex of halfedge _heh into its to-vertex. 
	 		ps->removeParticle(nvh.idx());
		}
	}
}

void SkinMeshShape::edgeFlip(ParticleOpenMesh &mesh)
{
	// clear status bits
	for(ParticleOpenMesh::EdgeIter e_it=mesh.edges_begin(); e_it!=mesh.edges_end();++e_it)
		e_it->tag = false;

	// swap to maximize minimum angle
	for(ParticleOpenMesh::EdgeIter e_it=mesh.edges_begin(); e_it!=mesh.edges_end();++e_it)
	{
		// do not flip on two neighboring edges
		
		if (e_it->tag)//changed by Matei Stroila
			continue;
		ParticleOpenMesh::EdgeHandle eh=e_it.handle();
		ParticleOpenMesh::HalfedgeHandle heRight=mesh.halfedge_handle(eh,0);
		ParticleOpenMesh::HalfedgeHandle heLeft=mesh.halfedge_handle(eh,1);
		ParticleOpenMesh::VertexHandle vhLeft=mesh.from_vertex_handle(heRight);
		ParticleOpenMesh::VertexHandle vhRight=mesh.from_vertex_handle(heLeft);
		ParticleOpenMesh::HalfedgeHandle heTopRight=mesh.next_halfedge_handle(heRight);
		ParticleOpenMesh::HalfedgeHandle heTopLeft=mesh.next_halfedge_handle(heTopRight);
		ParticleOpenMesh::HalfedgeHandle heBottomLeft=mesh.next_halfedge_handle(heLeft);
		ParticleOpenMesh::HalfedgeHandle heBottomRight=mesh.next_halfedge_handle(heBottomLeft);
		ParticleOpenMesh::VertexHandle vhTop=mesh.to_vertex_handle(heTopRight);
		ParticleOpenMesh::VertexHandle vhBottom=mesh.to_vertex_handle(heBottomLeft);
		
		double oldMinAngle=findMinAngle(mesh.point(vhTop),mesh.point(vhLeft),mesh.point(vhRight),mesh.point(vhBottom));
		double newMinAngle=findMinAngle(mesh.point(vhLeft),mesh.point(vhBottom),mesh.point(vhTop),mesh.point(vhRight));
		if (oldMinAngle*1.5<newMinAngle)
		{
			// flip
			if (mesh.is_flip_ok(eh))
			{
				mesh.flip(eh);
				// mark neighbors tags
				mesh.edge(eh).tag =true;
				mesh.edge(mesh.edge_handle(heTopRight)).tag = true;
				mesh.edge(mesh.edge_handle(heTopLeft)).tag = true;
				mesh.edge(mesh.edge_handle(heBottomLeft)).tag = true;
				mesh.edge(mesh.edge_handle(heBottomRight)).tag = true;
			}
		}
	}
}

/// triangle edge flips, splits or collapses
void SkinMeshShape::cleanup()
{
	// position must be a ParticleMesh.
	ParticleMesh *pm=dynamic_cast<ParticleMesh *>(position);
	if (!pm)
		return;
	ParticleOpenMesh &mesh=pm->mesh;

	// edge flips
	if (allowFlip)
		edgeFlip(mesh);
	edgeCollapseAndSplit(mesh);
}

double SkinMeshShape::findMinAngle(const ParticleOpenMesh::Point &a, const ParticleOpenMesh::Point &b,
		const ParticleOpenMesh::Point &c, const ParticleOpenMesh::Point &d)
{
	double a1=findMinAngle(a,b,c);
	double a2=findMinAngle(b,c,d);
	return (a1<a2?a1:a2);
}

double SkinMeshShape::findMinAngle(const ParticleOpenMesh::Point &a, const ParticleOpenMesh::Point &b,
		const ParticleOpenMesh::Point &c)
{
	double angleA=angleOf(b-a,c-a);
	double angleB=angleOf(c-b,a-b);
	double angleC=angleOf(a-c,b-c);

	double angle;
	if (angleA<angleB && angleA<angleC)
		angle=angleA;
	else if (angleB<angleA && angleB<angleC)
		angle=angleB;
	else
		angle=angleC;
	return angle;
}

double SkinMeshShape::angleOf(const ParticleOpenMesh::Point &a, const ParticleOpenMesh::Point &b)
{
	return std::acos((a[0]*b[0]+a[1]*b[1]+a[2]*b[2])/(a.length()*b.length()));
}

/** the new position is a weighted sum of
0.4*current value+0.4*neighbor length+0.2*target length
*/
void SkinMeshShape::calculateTargetEdgeLength()
{
		// position must be a ParticleMesh.
	ParticleMesh *pm=dynamic_cast<ParticleMesh *>(position);
	if (!pm)
		return;
	ParticleOpenMesh &mesh=pm->mesh;

	// for each particle
	for(ParticleOpenMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end();++v_it)
	{
		ParticleOpenMesh::VertexHandle vh=v_it.handle();
		ParticleOpenMesh::Point currentPoint(mesh.point(vh));
		// find neighbor length
		double neighborLength=0.0;
		int valence=0;
		for(ParticleOpenMesh::VertexVertexIter vv_it=mesh.vv_iter(vh); vv_it; ++vv_it)
		{
			ParticleOpenMesh::Point edge(currentPoint-mesh.point(vv_it.handle()));
			neighborLength+=edge.length();
			++valence;
		}
		neighborLength=neighborLength/valence;
	}

	// set new position
}
