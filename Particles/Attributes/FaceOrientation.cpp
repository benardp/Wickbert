/*
@file FaceOrientation.cpp
@author John C. Hart
@date 14 Apr. 2005
*/

#include "FaceOrientation.h"
#include "ParticleNormal.h"
#include "MeshInterrogator.h"
#include "Surface/OpenMesh/SurfaceMesh.h"
#include "pstools.h"
#include <time.h>

REGISTER_PARTICLESTUFF(FaceOrientation,"Attribute:FaceOrientation");

FaceOrientation::FaceOrientation(Particles *ps, const std::string& name)
	: ParticleNormal(ps,name)
{
	new Attached<MeshInterrogator>(this,&mi);
	new PSParamInt(this, &_maxNbParticles, 0,"nbParticles","mx nbParticles","this value indicates the maximum number of particles that will get an orientation on the faces in a random manor. 0 indicates a particle for all faces. Realize, that it remains compatible to FaceOrientation, as it uses the same randomization! Only the same amount of particles should be chosen.");
}

void FaceOrientation::attachAttributes()
{
	ParticleNormal::attachAttributes();

	if (!mi) return;

	SurfaceMesh *mesh = mi->getMesh();
	if (!mesh) return;

	mesh->update_face_normals();




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
			if (i >= ps->size()) {
				ps->addParticle();
			}
			SurfaceMesh::Normal nm(mesh->normal(faceHandles[i]));
			
			gmVector3 n(nm[0], nm[1], nm[2]);
			
			setNormal(i, n);
		}	
	}
	else
	{
		unsigned int i = 0;
		for (SurfaceMesh::FaceIter f_it = mesh->faces_begin(); f_it!=mesh->faces_end(); ++f_it) {
			if (i >= ps->size()) {
				ps->addParticle();
			}
			SurfaceMesh::Normal nm(mesh->normal(f_it.handle()));
			
			gmVector3 n(nm[0], nm[1], nm[2]);
			
			setNormal(i, n);
			i++;
		}
	}
}
