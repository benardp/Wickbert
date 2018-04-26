/*
@file ClusterMeshInterrogator.cpp
@author Wen Su
*/

#include "ClusterMeshInterrogator.h"
#include "ParticlePosition.h"
#include <ctime>
#include <limits>
#include <queue>
#include <functional>

REGISTER_PARTICLESTUFF(ClusterMeshInterrogator,"Behavior:ClusterMeshInterrogator");

ClusterMeshInterrogator::ClusterMeshInterrogator(Particles *ps)
	: ParticleBehavior(ps, std::string("ClusterMeshInterrogator"))
{
	iteration=10;
	minClusters=2;
	maxClusters=5;
	restartAlgorithm=0;
	phase=INITIALIZE;
}

void ClusterMeshInterrogator::attachAttributes()
{
	// find the cluster mash
	attachAttribute(surInt,std::string("SurfaceInterrogator"));

	ParticleBehavior::attachAttributes();
}

void ClusterMeshInterrogator::cleanup()
{
	// number of iteration for the clustering algorithm
	if (iteration==0)
		return;

	if (position==NULL)
		return;
	if (surInt==NULL)
		return;
	ClusterMesh *cm=dynamic_cast<ClusterMesh *>(surInt->surface);
	if (cm==NULL)
		return;
	ClusterOpenMesh &mesh=cm->mesh;

	if (restartAlgorithm)
	{
		restartAlgorithm=0;
		ps->clear();
		clusterCenters.clear();
		phase=INITIALIZE;
	}

	// four phases of the algorithm
	if (phase==INITIALIZE)
	{
		// assign a default center randomly
		srand(time(NULL));
		for(unsigned int i=0;i<minClusters;i++)
		{
			unsigned int j=rand()%mesh.n_faces();
			ClusterOpenMesh::FaceHandle fh(j);
			clusterCenters.push_back(fh);
		}
		// next phase
		phase=ASSIGNNEWCENTERS;
	}
	else if (phase==ASSIGNNEWCENTERS)
	{
		cm->removeCluster();
		lastFaceInCluster.clear();
		costHeap.clear();
		// clear the old positions
		ps->clear();
		// reinsert the new clusters
		for(unsigned int i=0;i<clusterCenters.size(); i++)
		{
			ClusterOpenMesh::FaceHandle &fh=clusterCenters[i];
			ClusterOpenMesh::Point &c(mesh.face(fh).center);
			unsigned int j=ps->addParticle();
			position->setPosition(j, gmVector3(c[0],c[1],c[2]));
			// change in heap
			costHeap.push_back(fh);
			lastFaceInCluster.push_back(fh);
			// assign cluster
			mesh.face(fh).cost=0;
			mesh.face(fh).cluster=j;
		}
		// make heap really does not do anything because all cost is 0
		std::make_heap(costHeap.begin(),costHeap.end(),FaceCostLessThan<ClusterOpenMesh>(mesh));
		phase=EXPANDCLUSTERS;
	}
	else if (phase==EXPANDCLUSTERS)
	{
		if (costHeap.size()==0)
		{
			// this phase is done
			phase=FINDCLUSTERCENTERS;
		}
		else 
		{
			// fh is the current face
			std::pop_heap(costHeap.begin(),costHeap.end(),FaceCostLessThan<ClusterOpenMesh>(mesh));
			ClusterOpenMesh::FaceHandle fh=costHeap.back();
			lastFace=fh;
			lastFaceInCluster[mesh.face(fh).cluster]=fh;
			costHeap.pop_back();
			float currentCost=mesh.face(fh).cost;
			// for each edge neighbor of face
			ClusterOpenMesh::HalfedgeHandle endhe=mesh.halfedge_handle(fh);
			ClusterOpenMesh::HalfedgeHandle he=endhe;
			do
			{
				// if it is not boundary
				if (mesh.is_boundary(he))
					continue;
				// the neighbor face
				ClusterOpenMesh::Face &face=mesh.face(mesh.face_handle(mesh.opposite_halfedge_handle(he)));
				float newCost=currentCost+mesh.edge(mesh.edge_handle(he)).cost;
				if (newCost<face.cost)
				{
					face.cost=newCost;
					face.cluster=mesh.face(fh).cluster;
					costHeap.push_back(mesh.handle(face));
					std::push_heap(costHeap.begin(),costHeap.end(),FaceCostLessThan<ClusterOpenMesh>(mesh));
				}
				he=mesh.next_halfedge_handle(he);
			}
			while (he!=endhe);
		}
	}
	else if (phase==FINDCLUSTERCENTERS)
	{
		// next phase
		phase=ASSIGNNEWCENTERS;
		--iteration;
		// store new positions
		std::vector <ClusterOpenMesh::Point> center;
		std::vector <unsigned int> count;
		for(unsigned int i=0;i<ps->size();i++)
		{
			center.push_back(ClusterOpenMesh::Point(0,0,0));
			count.push_back(0);
		}
		// find new center
		for (ClusterOpenMesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end();++f_it)
		{
			int id=f_it->cluster;
			count[id]=count[id]+1;
			center[id]=center[id]+f_it->center;
		}
		// find average position
		for(unsigned int i=0;i<ps->size();i++)
			center[i]=center[i]/count[i];
		// find the new face center by finding the geomtric center of the faces
		for(unsigned int i=0;i<ps->size();i++)
		{
			float min=std::numeric_limits<float>::max();
			// all face and this center
			for (ClusterOpenMesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end();++f_it)
			{
				ClusterOpenMesh::Point vec(f_it->center-center[i]);
				float dist=vec.length();
				if (dist<min)
				{
					min=dist;
					clusterCenters[i]=f_it.handle();
				}
			}
		}
		// add the last face found as the new cluster center
		if (clusterCenters.size()<maxClusters)
		{
			clusterCenters.push_back(lastFace);
		}
		// print some debug info
		std::cout<<"iteration: "<<iteration<<std::endl;
		std::cout<<"last face costs: "<<mesh.face(lastFace).cost<<std::endl;
		for(unsigned int i=0;i<lastFaceInCluster.size();i++)
		{
			std::cout<<"cluster: "<<i<< " last face cost: " << mesh.face(lastFaceInCluster[i]).cost<<std::endl;
		}
	}
}

int ClusterMeshInterrogator::qlen()
{
	return 4;
}
void ClusterMeshInterrogator::getq(double *q)
{
	q[0] = minClusters;
	q[1] = maxClusters;
	q[2] = restartAlgorithm;
	q[3] = iteration;
}
void ClusterMeshInterrogator::setq(double *q)
{
	minClusters = q[0];
	maxClusters = q[1];
	restartAlgorithm = q[2];
	iteration = q[3];
}
void ClusterMeshInterrogator::qname(char **qn)
{
	qn[0] = "minClusters";
	qn[1] = "maxClusters";
	qn[2] = "restartAlgorithm";
	qn[3] = "iteration";
}
