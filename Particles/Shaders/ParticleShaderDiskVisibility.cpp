#include "ParticleShaderDiskVisibility.h"

REGISTER_PARTICLESTUFF(ParticleShaderDiskVisibility,"Shader:ParticleShaderDiskVisibility");

ParticleShaderDiskVisibility::ParticleShaderDiskVisibility(Particles *ps)
	:ParticleShaderDisk(ps)
{
	setName(std::string("ParticleShaderDiskVisibility"));

	new PSParamBool(this,&useVisibility,false,"useVisibility","Use Visibility","Turn it on if you want to use Visibility.");
	new Attached<ParticleVisibility>(this,&vis);
	new Attached<ViewDependence>(this,&view);
	
	cameraPosChanged = false;
	oldCameraPosition = new gmVector3(0.0,0.0,0.0);
}

void ParticleShaderDiskVisibility::attachAttributes()
{
	ParticleShaderDisk::attachAttributes();
}

void ParticleShaderDiskVisibility::drawPre()
{
	myCameraPosition = view->getCameraPosition();

	if( fabs((*oldCameraPosition)[0] - (*myCameraPosition)[0]) < CAMERAPOS_SENS
		&&
		fabs((*oldCameraPosition)[1] - (*myCameraPosition)[1]) < CAMERAPOS_SENS
		&&
		fabs((*oldCameraPosition)[2] - (*myCameraPosition)[2]) < CAMERAPOS_SENS
		)
	{
		cameraPosChanged = false;
		return;
	}
	else
	{
		*oldCameraPosition = *myCameraPosition;
		cameraPosChanged = true;
		return;
	}
}


void ParticleShaderDiskVisibility::drawParticle(int i)
{
	int isVisible = 1;
	if(useVisibility)
	{
		if(cameraPosChanged) 
			isVisible = vis->getVisibility(myCameraPosition,i);
		else
			isVisible = vis->getVisibility(i);
	}				
	if(isVisible)
		ParticleShaderDisk::drawParticle(i);
}
