//Only compiled if CG Shaders are enabled
#ifdef WB_USE_CG

/**
 * @file ParticleShaderVoronoi.cpp
 * @author Jerry O. Talton III
 */

#include "ParticleShaderVoronoiDevelop.h"
#include "AdaptiveRepulsionData.h"
#include "ParticleOrientation.h"
#include "ParticlePosition.h"

REGISTER_PARTICLESTUFF(ParticleShaderVoronoiDevelop,"Shader:ParticleShaderVoronoiDevelop");

ParticleShaderVoronoiDevelop::ParticleShaderVoronoiDevelop(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderVoronoiDevelop"))
{

	quad=gluNewQuadric();

	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleOrientation>(this,&orientation);

	new PSParamString(this,&radius_source,"AdaptiveRepulsionData:radius","radattr","radius source",
		"Source of the radius of the form Attribute:param");
	new PSParamDouble(this,&scale,1.35,"scale","scale factor",
		"Amount to scale the radius");
	new PSParamInt(this,&sides,10,"sides","polygon sides",
		"# of sides of the disk approximating polygons");
	new PSParamDouble(this,&toleranceScale,1.8,"tolerance scale",
		"Number of radii we limit fragments to");

	radius_attr = NULL;
	radius_data = NULL;

	pbufferOne = NULL;

	// Get the current viewport size
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	// Check to make sure the necessary extensions are supported
	pbufferOne->testExtensions();

	// Setup the pbuffers
	setupBuffers(viewport[2], viewport[3]);

	// Setup textures
	glGenTextures(1, &depth_texture);		
	glBindTexture(GL_TEXTURE_RECTANGLE_NV, depth_texture);
	glTexParameteri(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	// Setup the CG environment
	setupCG();

}

void ParticleShaderVoronoiDevelop::attachAttributes()
{

	PSParamPerParticle *radius_param = NULL;

	ParticleShader::attachAttributes();

	std::string attr_name = radius_source.substr(0,radius_source.find(':'));
	std::string param_name = radius_source.substr(radius_source.find(':')+1);
	if (radius_attr = ps->getAttribute<ParticleAttribute>(attr_name))
		radius_param = radius_attr->perparticle.findparam(param_name);
	if (radius_param)
		radius_data = (DoubleVector *)radius_param->ref();
	else
		radius_attr = NULL;
}

void ParticleShaderVoronoiDevelop::setupBuffers(int width, int height)
{
	hdc = wglGetCurrentDC();
	hglrc = wglGetCurrentContext();

	pbufferOne = new nv_pbuffer(width, height, 1);

	pbufferOne->wglGetLastError();	

	// If we're sharing contexts, don't need to share lists
	if (pbufferOne->onecontext == false) 
		wglShareLists(pbufferOne->hglrc, hglrc);
	
}

void ParticleShaderVoronoiDevelop::bindCGParametersVertex()
{

	// Bind the vertex program parameters
	modelViewProj = cgGetNamedParameter(vProgram,"modelViewProj");
	modelView = cgGetNamedParameter(vProgram,"modelView");
	modelViewIT = cgGetNamedParameter(vProgram,"modelViewIT");
	viewMatrix = cgGetNamedParameter(vProgram,"viewMatrix");
	diskCenterWC = cgGetNamedParameter(vProgram,"diskCenterWC");

	if (!modelViewProj || !modelView || !modelViewIT || !viewMatrix || !diskCenterWC)
		exit(-1);
}


void ParticleShaderVoronoiDevelop::bindCGParametersFragmentOne()
{

	// Bind the first fragment program parameters
	lightPositionEC = cgGetNamedParameter(fProgramOne,"lightPositionEC");
	diffuseColor = cgGetNamedParameter(fProgramOne,"diffuseColor");
	diskRadius = cgGetNamedParameter(fProgramOne,"diskRadius");
	firstPassDepth = cgGetNamedParameter(fProgramOne,"firstPassDepth");
	depthRange = cgGetNamedParameter(fProgramOne,"depthRange");	
	toleranceScaleFactor = cgGetNamedParameter(fProgramOne,"toleranceScaleFactor");

	cgGLSetTextureParameter(firstPassDepth, depth_texture);

	if (!lightPositionEC || !diffuseColor || !diskRadius || !firstPassDepth
		|| !depthRange || !toleranceScaleFactor)
		exit(-1);
}

void ParticleShaderVoronoiDevelop::setupCG()
{

	if (cgGLIsProfileSupported(CG_PROFILE_VP30))
		vProfile = CG_PROFILE_VP30;
	else exit(-1);

	if (cgGLIsProfileSupported(CG_PROFILE_FP30))
		fProfile = CG_PROFILE_FP30;
	else exit(-1);

	// Create the Cg context
	CGcontext cgContext = cgCreateContext();
	
	// Create the vertex program
	vProgram = cgCreateProgramFromFile(cgContext, CG_SOURCE, "../Particles/vertex.cg", 
		vProfile, NULL, NULL);
	cgGLLoadProgram(vProgram);

	// Create the fragment programs
	fProgramOne = cgCreateProgramFromFile(cgContext, CG_SOURCE, "../Particles/accum_frag_d.cg",
		fProfile, NULL, NULL);
	cgGLLoadProgram(fProgramOne);

}

void ParticleShaderVoronoiDevelop::drawPost()
{


	// Get the viewport dimensions
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	if ((pbufferOne->height != viewport[3]) || (pbufferOne->width != viewport[2]))
	{
		delete pbufferOne;
		setupBuffers(viewport[2], viewport[3]);
	}

	// Get the light position
	GLfloat lightpos[4];
	glGetLightfv(GL_LIGHT0, GL_POSITION, lightpos);

	// Create storage for matrices
	GLfloat mm[16];
	GLfloat pm[16];

	// If we're not sharing contexts, need to copy state
	glGetFloatv(GL_MODELVIEW_MATRIX, mm);
	glGetFloatv(GL_PROJECTION_MATRIX, pm);

	// Switch to pbufferOne rendering context
	if (wglMakeCurrent(pbufferOne->hdc, pbufferOne->hglrc) == FALSE)
		pbufferOne->wglGetLastError();

	// If we're not sharing contexts, need to update state
	if (pbufferOne->onecontext == false) {

		glMatrixMode(GL_PROJECTION);
		glLoadMatrixf(pm);
	
		glMatrixMode(GL_MODELVIEW);
		glLoadMatrixf(mm);

		glEnable(GL_DEPTH_TEST);
	}
	else {
		// Turn off color writes, for better performance
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	}

	// Clear the depth texture
	glClear(GL_DEPTH_BUFFER_BIT);

	// Draw into the buffer for the depth texture
	for (unsigned int i = 0; i < ps->size(); i++)
		drawParticle(i);

	// Turn back on color for eventual rendering
	if (pbufferOne->onecontext) {
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	}

	// Switch to normal rendering context
	if (wglMakeCurrent(hdc, hglrc) == FALSE)
		pbufferOne->wglGetLastError();
	
	// Setup CG
	cgGLEnableProfile(vProfile);
	cgGLEnableProfile(fProfile);
	cgGLBindProgram(vProgram);
	bindCGParametersVertex();
	bindCGParametersFragmentOne();
	cgGLBindProgram(fProgramOne);
	
	// Make the depth texture active
	glActiveTextureARB = (PFNGLACTIVETEXTUREARBPROC) wglGetProcAddress("glActiveTextureARB");
	glActiveTextureARB(GL_TEXTURE0_ARB);
	
	// Bind the depth texture
	glEnable(GL_TEXTURE_RECTANGLE_NV);
	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glBindTexture(GL_TEXTURE_RECTANGLE_NV, depth_texture);
	
	// Bind pbufferOne as a depth texture
	wglBindTexImageARB = (PFNWGLBINDTEXIMAGEARBPROC)wglGetProcAddress("wglBindTexImageARB");
	if (wglBindTexImageARB(pbufferOne->hpbuffer, WGL_DEPTH_COMPONENT_NV) == FALSE)
		pbufferOne->wglGetLastError();

	cgGLEnableTextureParameter(firstPassDepth);
	
	glActiveTextureARB(GL_TEXTURE1_ARB);
	
	// Set parameters
	cgGLSetParameter1f(depthRange, GLfloat(abs(10.0 - 0.5)));
	cgGLSetParameter4fv(lightPositionEC, lightpos);
	cgGLSetParameter1f(toleranceScaleFactor, GLfloat(toleranceScale));
	cgGLSetStateMatrixParameter(viewMatrix, CG_GL_MODELVIEW_MATRIX, CG_GL_MATRIX_IDENTITY);

	// Clear the normal texture
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
	
	// Turn on blending for accumulation
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ONE);
	glDisable(GL_DEPTH_TEST);

	// Increase the scale of disks so as to create more overlap
	double tempscale = scale;
	scale = 2.0 * tempscale;

	// Draw the particles to the screen
	for (unsigned int i = 0; i < ps->size(); i++) {
		drawSecondPass(i);
	}

	// Reset scale
	scale = tempscale;

	glEnable(GL_DEPTH_TEST);
	// Disable blending
	glDisable(GL_BLEND);

	cgGLDisableTextureParameter(firstPassDepth);
	
	glDisable(GL_TEXTURE_RECTANGLE_NV);

	// Release the depth texture
	wglReleaseTexImageARB = (PFNWGLRELEASETEXIMAGEARBPROC)wglGetProcAddress("wglReleaseTexImageARB");
	if (wglReleaseTexImageARB(pbufferOne->hpbuffer, WGL_DEPTH_COMPONENT_NV) == false)
		pbufferOne->wglGetLastError();

	cgGLDisableProfile(vProfile);
    cgGLDisableProfile(fProfile);
}



void ParticleShaderVoronoiDevelop::drawParticle(int i)
{
	// push the shader name down
	glPushName(0xffffffff);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	gmVector3 p = position->getPosition(i);
	glTranslatef(p[0],p[1],p[2]);

	gmMatrix4 rotMat = orientation->getMatrix(i);

	GLdouble mat[16];
	rotMat.copyTo(mat);

	glMultMatrixd(mat);	

	drawShape(i);
	
	glPopMatrix();
	glPopName();
}

void ParticleShaderVoronoiDevelop::drawSecondPass(int i)
{

	// push the shader name down
	glPushName(0xffffffff);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	gmVector3 p = position->getPosition(i);
	glTranslatef(p[0],p[1],p[2]);

	gmMatrix4 rotMat = orientation->getMatrix(i);

	GLdouble mat[16];
	rotMat.copyTo(mat);

	glMultMatrixd(mat);

	//setColor(i);

	// Pull the diffuse color
	GLfloat diff[4];
	glGetMaterialfv(GL_FRONT, GL_DIFFUSE, diff);

	// Set the disk center
	GLfloat dC[3] = {0.0, 0.0, 0.0};

	// Get the radius
	if (radius_data)
		radius = (*radius_data)[i];
	else
		radius = 1.0;
	
	// Set the per-particle Cg parameters
	cgGLSetParameter3fv(diskCenterWC, dC);
	cgGLSetParameter4fv(diffuseColor, diff);
	cgGLSetParameter1f(diskRadius, GLfloat(radius * scale));
	cgGLSetStateMatrixParameter(modelViewProj, CG_GL_MODELVIEW_PROJECTION_MATRIX, CG_GL_MATRIX_IDENTITY);
	cgGLSetStateMatrixParameter(modelView, CG_GL_MODELVIEW_MATRIX, CG_GL_MATRIX_IDENTITY);
	cgGLSetStateMatrixParameter(modelViewIT, CG_GL_MODELVIEW_MATRIX, CG_GL_MATRIX_INVERSE_TRANSPOSE);
	
	// Now, actually draw the shape
	drawShape(i);	
   
	glPopMatrix();
	glPopName();

}


void ParticleShaderVoronoiDevelop::drawShape(int i)
{

	if (radius_data)
		radius = (*radius_data)[i];
	else
		radius = 1.0;

	// Turn on backface culling, since it's integral to the algorithm
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	gluDisk(quad,0,radius*scale,sides,sides);

	glDisable(GL_CULL_FACE);
}


ParticleShaderVoronoiDevelop::~ParticleShaderVoronoiDevelop()
{
	if (pbufferOne) delete pbufferOne;
	gluDeleteQuadric(quad);
	cgDestroyContext(cgContext);
}

#endif