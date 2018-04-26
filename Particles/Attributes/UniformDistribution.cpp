/**
 * @file ParticleVelocity.cpp
 * @author Elmar Eisemann
 */

#include "UniformDistribution.h"
#include "ParticlePosition.h"
#include "ParticleBoundingBox.h"

REGISTER_PARTICLESTUFF(UniformDistribution,"Attribute:UniformDistribution");

UniformDistribution::UniformDistribution(Particles *ps, const std::string& name)
	: ParticleAttribute(ps, name) 
{
	new Attached<ParticlePosition>(this, &_thisPosition);	

	new PSParamInt(this, &(_nbParticles[0]), 40,"nbX", "particles along X", "Number of particles along the x axis.");
	new PSParamInt(this, &(_nbParticles[1]), 40,"nbY", "particles along Y", "Number of particles along the y axis.");
	new PSParamInt(this, &(_nbParticles[2]), 40,"nbZ", "particles along Z", "Number of particles along the z axis.");

	new PSAttrRefParam(this, &_boundingBoxAttributeRef, "invalid", "bbox", "bounding box", "Assign the bounding box of a particle system to bound distribution");

	new PSParamBool(this, &_sampleMaxBorder,true,"maxB","include Max", "Include the maximum border as a sampling point, this can not be wanted in the case of nearest neighbor filtering");

	new PSParamButton(this, new ApplyUniformRedistributionCallback(this),"redistrib.","redistribute",
		"Distribute particles in bounding box.");
}

void UniformDistribution::redistributeUniformly()
{
	ParticleBoundingBox* boundingBox=dynamic_cast<ParticleBoundingBox*>(_boundingBoxAttributeRef);
	if (boundingBox==0)
		return;
	const gmVector3 & minC=boundingBox->min;
	const gmVector3 & maxC=boundingBox->max;

	const unsigned int nbParticles = ps->size();
	const unsigned int nbWanted=_nbParticles[0]*_nbParticles[1]*_nbParticles[2];
	if (nbParticles<nbWanted)
	{
		for (unsigned int i=0;i<nbWanted-nbParticles;++i)
		{
			ps->addParticle();
		}
	}
	else if (nbParticles>nbWanted)
	{
		//remove all the last particles, until the size is ok
		//this is the way it is supposed to work... although 
		//I do personnally not like the fact that when you delete 
		//a particle in the middle, the indices change
		for (unsigned int i=0;i<nbParticles-nbWanted;++i) 
			ps->removeParticle(nbParticles-i-1); 
	}
	assert(ps->size()==nbWanted);


	gmVector3 offset;

	if (_sampleMaxBorder)
	{
		offset = gmVector3(	(maxC[0]-minC[0])/(_nbParticles[0]-1),
							(maxC[1]-minC[1])/(_nbParticles[1]-1),
							(maxC[2]-minC[2])/(_nbParticles[2]-1)
						);
	}
	else
	{
		offset = gmVector3(	(maxC[0]-minC[0])/(_nbParticles[0]),
							(maxC[1]-minC[1])/(_nbParticles[1]),
							(maxC[2]-minC[2])/(_nbParticles[2])
						);
	}

	for (unsigned int z=0;z<(unsigned int) _nbParticles[2];++z)
		for (unsigned int y=0;y<(unsigned int) _nbParticles[1];++y)
			for (unsigned int x=0;x<(unsigned int) _nbParticles[0];++x)
			{	
				_thisPosition->setPosition(
											x+_nbParticles[0]*y+_nbParticles[0]*_nbParticles[1]*z,
											minC+gmVector3(offset[0]*x,offset[1]*y,offset[2]*z)
											);
			}
}
