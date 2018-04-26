/*
@file FacePosition.cpp
@author John C. Hart
@date 4 Jan. 2005
*/

#include "FacePosition.h"
#include "ParticlePosition.h"
#include "MeshInterrogator.h"
#include "Surface/OpenMesh/SurfaceMesh.h"
#include "pstools.h"
#include <time.h>

REGISTER_PARTICLESTUFF(FacePosition,"Attribute:FacePosition");

FacePosition::FacePosition(Particles *ps, const std::string& name)
	: ParticlePosition(ps,name)
{
	// new Attached<ParticlePosition>(this,&position);
	new Attached<MeshInterrogator>(this,&mi);

	new PSParamInt(this, &_maxNbParticles, 0,"nbParticles","mx nbParticles","this value indicates the maximum number of particles that will be placed on the faces in a random manor. 0 indicates on all faces. Realize, that it remains compatible to FaceOrientation, as it uses the same randomization! Only the same amount of particles should be chosen.");
}

void FacePosition::attachAttributes()
{
	ParticlePosition::attachAttributes();

	if (!mi) return;

	SurfaceMesh *mesh = mi->getMesh();
	if (!mesh) return;


	if (_maxNbParticles>0)
	{
		std::vector<FaceHandle> faceHandles;
		faceHandles.reserve(mesh->n_faces());
		for (SurfaceMesh::FaceIter f_it = mesh->faces_begin(); f_it!=mesh->faces_end(); ++f_it) 
		{
			faceHandles.push_back(f_it.handle());	
		}

		//we always want the same random sequence.
		srand(STANDARD_RANDOM_SEED);
		std::pointer_to_unary_function<int,int> pt(posIntRand);
		std::random_shuffle(faceHandles.begin(),faceHandles.end(), pt);
		srand(time(NULL));
		const unsigned int minimum=(unsigned int)(((unsigned int)_maxNbParticles)<mesh->n_faces())?_maxNbParticles:mesh->n_faces();

		for (unsigned int i=0; i<minimum; ++i)
		{

			unsigned int n=0;
			SurfaceMesh::Point bary(0,0,0);
			for (SurfaceMesh::FaceVertexIter fv_it = mesh->fv_iter(faceHandles[i]); fv_it; ++fv_it) {
				SurfaceMesh::Point v(mesh->point(fv_it));
				bary += v;
				n++;
			}
			bary *= (float) 1.0/n;
			if (i >= ps->size()) 
				ps->addParticle();
			x[i][0] = bary[0];
			x[i][1] = bary[1];
			x[i][2] = bary[2];
		}	
	}
	else
	{

		unsigned int i = 0;
		for (SurfaceMesh::FaceIter f_it = mesh->faces_begin(); f_it!=mesh->faces_end(); ++f_it) {
			int n = 0;
			SurfaceMesh::Point bary(0,0,0);
			for (SurfaceMesh::FaceVertexIter fv_it = mesh->fv_iter(f_it.handle()); fv_it; ++fv_it) {
				SurfaceMesh::Point v(mesh->point(fv_it));
				bary += v;
				n++;
			}
			bary *= (float) 1.0/n;
			if (i >= ps->size()) 
				ps->addParticle();
			x[i][0] = bary[0];
			x[i][1] = bary[1];
			x[i][2] = bary[2];

			i++;
		}
	}
}
