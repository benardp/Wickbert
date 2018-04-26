/*
@file ParticleShaderDisk.h
@author Wen Su
This class draws a disk for a particle object.
*/

#ifndef PARTICLESHADERDISKDEVELOP_H
#define PARTICLESHADERDISKDEVELOP_H

#include "ParticleShader.h"
#include <CG/cg.h>
#include <CG/cgGL.h>

class ParticleShaderDiskDevelop: public ParticleShader
{

protected:
	GLUquadricObj* quad;
	double radius;
	double scale;
	int sides;

	std::string radius_source;
	ParticleAttribute *radius_attr;
	DoubleVector *radius_data;

	// this gives informatio for shader to draw the size of objects
	AdaptiveRepulsionData *ARData;

	CGcontext cgContext;				// The Cg context
	CGprogram vProgram;					// The vertex program
	CGprogram fProgram;					// The fragment programs
	CGprofile vProfile, fProfile;		// The vertex and fragment profiles

	CGparameter modelView, modelViewIT;	// The modelview matrix
	CGparameter modelViewProj;			// The modelview matrix and the projection matrix
	CGparameter viewMatrix;				// A matrix containing just the camera transform

	CGparameter diffuseColor;			// The diffuse color of the particles
	CGparameter lightPositionEC;		// The position of the light in eye coordinates
	
	CGparameter diskCenterWC;			// The center of the each disk in world coordinate
	CGparameter diskRadius;				// The radius of the disks

public:

	MAKE_PARTICLESTUFF_NAME();

	ParticleShaderDiskDevelop(Particles *ps=NULL);
	virtual void attachAttributes();
   
	virtual ~ParticleShaderDiskDevelop();

	void setupCG();

	void bindCGParametersFragment();
	void bindCGParametersVertex();

	virtual void drawParticle(int i);
};

#endif