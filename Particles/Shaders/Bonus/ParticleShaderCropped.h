/*
@file ParticleShaderCropped.h
@author Jerry O. Talton II
This class draws a cropped disk for each particle object.
*/

#ifndef PARTICLESHADERCROPPED_H
#define PARTICLESHADERCROPPED_H

#include "ParticleShader.h"
#include "ImplicitInterrogator.h"
#include "nv_pbuffer.h"
#include <CG/cg.h>
#include <CG/cgGL.h>


class ParticleShaderCropped: public ParticleShader
{

protected:

	double radius;
	double scale;
	int sides;

	std::string radius_source;
	ParticleAttribute *radius_attr;
	DoubleVector *radius_data;

	double toleranceScale;

	GLUquadricObj* quad;

	/// this gives information for shader to draw the size of objects
	AdaptiveRepulsionData *ARData;

	nv_pbuffer *pbuffer;				// The pbuffer, used for the depth pass
	
	HDC hdc;							// The device context
	HGLRC hglrc;						// The OpenGL rendering context

	GLuint depth_texture;				// Handle to our depth texture

	/* Cg variables for passing data to our vertex and fragment programs */

	CGcontext cgContext;				// The Cg context
	CGprogram vProgram, fProgram;		// The vertex and fragment programs
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

	ImplicitInterrogator *imp_int;


public:


	MAKE_PARTICLESTUFF_NAME();

	ParticleShaderCropped(Particles *ps=NULL);
	
	virtual void attachAttributes();
   
	virtual ~ParticleShaderCropped();

	/// The function that sets up the pbuffer
	void setupBuffer(int height, int width);

	/// The function that sets up the Cg environment
	void setupCG();

	void bindCGParameters();

	virtual void draw(int s);

	virtual void drawParticle(int i, int s);

	/// An overloaded version of gluDisk(), to allow for greater customization
	/// and flexibility.
	virtual void drawShape(int i);

	/// Draws the second pass with cg enabled
	void drawSecondPass(int i, int s);

	/// parameters
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);

};

#endif
