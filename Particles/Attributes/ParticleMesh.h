/*
@file ParticleMesh.h
@author Wen Su
Particles are vertices in a triangle mesh
*/

#ifndef PARTICLEMESH_H
#define PARTICLEMESH_H

#include <string>
#include "ParticlePosition.h"
#include "ParticleAttribute.h"

#ifdef check
#undef check
#endif

#include "Surface/OpenMesh/ParticleMeshTraits.h"

/** ParticleMesh is a ParticleAttribute that stores a mesh data structure.
 * We currently use the OpenMesh library and the mesh is maintained in the
 * member variable mesh of type ParticleOpenMesh. This attribute inherits
 * ParticlePosition and the members get/setPosition maintain the position
 * of the mesh vertices.
 */
class ParticleMesh : public ParticlePosition
{
public:
	/// OpenMesh mesh data structure 
	ParticleOpenMesh mesh;

	/// Filename used to read in a mesh from a source file
	std::string filename;
    std::string outFileName;

private:
	/// Filename of mesh currently loaded. Does not reload if it matches filename.
	std::string loaded_filename;

    class WriteMeshCB : public PSParamButton::Callback {

    public:
        WriteMeshCB(ParticleMesh *mesh);
        virtual void onbuttonpress();
    private:
        ParticleMesh *_mesh;
    };

public:
	MAKE_PARTICLESTUFF_NAME();

	/// default constructor
	ParticleMesh(Particles *ps=NULL);

	/// create diamond
	void createDiamond(double size=10.0f);

	/// this gives a calling convertion for children to get the position of index i
	virtual gmVector3 getPosition(unsigned int i);
	
	/// this gives a calling convertion for children to set the position of index i
	virtual void setPosition(unsigned int i, gmVector3 p);

	virtual void clear();

	virtual void attachAttributes();

	/// Callback for particle addition.
	void virtual particleAdded();

	/// Callback for particle removal.
	void virtual particleRemoved(unsigned int i);

    void writeMesh();
};

#endif
