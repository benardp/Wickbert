#include "MeshInterrogator.h"
#include "Surface/Surface.h"
#include "Surface/OpenMesh/SurfaceMesh.h"
#include "ParticleSystem.h"

REGISTER_PARTICLESTUFF(MeshInterrogator,"Attribute:MeshInterrogator");

MeshInterrogator::MeshInterrogator(Particles* ps, SurfaceMesh* ms, const std::string& name)
	: SurfaceInterrogator(ps, ms, name)
{
}

void MeshInterrogator::attachAttributes()
{
	Surfaces *surfs = ps->particleSystem->surfaces;
	if (surfs)
		surface = dynamic_cast<SurfaceMesh *>(surfs->withName(surfname));
}

void MeshInterrogator::setSurface(Surface *s)
{
	surface = dynamic_cast<SurfaceMesh *>(s);
}

void MeshInterrogator::setMesh(SurfaceMesh *ms)
{
	surface = ms;
}

SurfaceMesh* MeshInterrogator::getMesh()
{
	return dynamic_cast<SurfaceMesh *>(surface);
}
