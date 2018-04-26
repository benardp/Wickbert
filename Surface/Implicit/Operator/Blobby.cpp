/** @file Convolution.cpp 
 * @author John C. Hart
 * @date 10 Jan. 2006
 */

#include "Blobby.h"

REGISTER_IMPLICIT(Blobby,"Blobby");

Blobby::Blobby(void)
{
	//MEMORY LEAK ??? -- EE
	//but this is probably not a real class, I guess...
	basis = new Quadric(1.0,1.0,1.0,-1.0);
	blinn = new Blinn(basis,1.0);

	new SurfParamDouble(this,&blobbiness,1.0,"B","blobbiness",
		"Blobbiness of blob blending.");

	new SurfAttrRefParam(this,(ParticleAttribute **)&positions,"ParticlePosition","pos","positions",
		"Centers of blobs.");

	new SurfAttrRefParam(this,(ParticleAttribute **)&radii,"ParticleScalar","rad","radii",
		"Radii of blobs.");
}

