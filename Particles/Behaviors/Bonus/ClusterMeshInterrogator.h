/*
@file ClusterMeshInterrogator.h
@author Wen Su
Particles are clusters on the mesh
*/

#ifndef CLUSTERMESHINTERROGATOR_H
#define CLUSTERMESHINTERROGATOR_H

#include "Particles.h"
#include "SurfaceInterrogator.h"
#include "ParticleBehavior.h"
#include "Surface/OpenMesh/ClusterMesh.h"
#include <vector>

class ClusterMeshInterrogator : public ParticleBehavior
{
public:
	MAKE_PARTICLESTUFF_NAME();

	/// different phase of the algorithm
	enum ClusterPhase {INITIALIZE, ASSIGNNEWCENTERS, EXPANDCLUSTERS, FINDCLUSTERCENTERS};
	ClusterPhase phase;

	/// the position is stored in particle position this stores the face handles.
	std::vector <ClusterOpenMesh::FaceHandle> clusterCenters;

	/// this vector will be used as a heap to do the Dijkstra walk
	std::vector <ClusterOpenMesh::FaceHandle> costHeap;

	/// stores the last face that is added to this cluster
	std::vector <ClusterOpenMesh::FaceHandle> lastFaceInCluster;

	/// last face that is added to the last iteration
	ClusterOpenMesh::FaceHandle lastFace;

	/// parameters
	unsigned int iteration;
	unsigned int minClusters;
	unsigned int maxClusters;
	int restartAlgorithm;

	/// pointer to the surface mesh
	SurfaceInterrogator* surInt;

	/// default constructor
	ClusterMeshInterrogator(Particles *ps=NULL);

	/// cleanup does the clustering step
	virtual void cleanup();

	virtual void attachAttributes();

	/// parameters
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

};

#endif
