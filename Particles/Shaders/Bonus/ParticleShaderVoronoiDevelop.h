/*
@file ParticleShaderVoronoi.h
@author Jerry O. Talton II
This class draws a cropped, Voronoi shaded disk for each particle object.
*/

#ifndef PARTICLESHADERVORONOIDEVELOP_H
#define PARTICLESHADERVORONOIDEVELOP_H

#include "ParticleShader.h"
#include "nv_pbuffer.h"
#include <CG/cg.h>
#include <CG/cgGL.h>


class ParticleOrientation;
class ParticlePosition;


class ParticleShaderVoronoiDevelop: public ParticleShader
{

protected:

	double radius;
	double scale;
	int sides;

	ParticlePosition *position;
	ParticleOrientation *orientation;

	std::string radius_source;
	ParticleAttribute *radius_attr;
	DoubleVector *radius_data;

	double toleranceScale;

	GLUquadricObj* quad;

	/// this gives information for shader to draw the size of objects
	AdaptiveRepulsionData *ARData;

	nv_pbuffer *pbufferOne;				// The pbuffer, used for the depth pass
	
	HDC hdc;							// The device context
	HGLRC hglrc;						// The OpenGL rendering context

	GLuint depth_texture;				// Handle to our depth texture

	/* Cg variables for passing data to our vertex and fragment programs */

	CGcontext cgContext;				// The Cg context
	CGprogram vProgram;					// The vertex program
	CGprogram fProgramOne, fProgramTwo; // The fragment programs
	CGprofile vProfile, fProfile;		// The vertex and fragment profiles

	CGparameter modelView, modelViewIT;	// The modelview matrix
	CGparameter modelViewProj;			// The modelview matrix and the projection matrix
	CGparameter viewMatrix;				// A matrix containing just the camera transform
	CGparameter depthRange;				// The length of the viewing frustrum
	
	CGparameter diffuseColor;			// The diffuse color of the particles

	CGparameter lightPositionEC;		// The position of the light in eye coordinates
	CGparameter diskCenterWC;			// The center of the each disk in world coordinate
	
	CGparameter diskRadius;				// The radius of the disks
	CGparameter toleranceScaleFactor;	// The amount to multiply the radius by for the tolerance
	
	CGparameter firstPassDepth;			// The depth texture


public:


	MAKE_PARTICLESTUFF_NAME();

	ParticleShaderVoronoiDevelop(Particles *ps=NULL);
	
	virtual void attachAttributes();
   
	virtual ~ParticleShaderVoronoiDevelop();

	/// The function that sets up the pbuffer
	void setupBuffers(int height, int width);

	/// The function that sets up the Cg environment
	void setupCG();

	void bindCGParametersVertex();
	void bindCGParametersFragmentOne();

	virtual void drawPost();

	virtual void drawParticle(int i);

	virtual void drawShape(int i);

	/// Draws the second pass with cg enabled
	void drawSecondPass(int i);

};

#endif

