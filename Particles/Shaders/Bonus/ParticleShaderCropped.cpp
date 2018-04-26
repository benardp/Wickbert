#include "ParticleShaderCropped.h"
#include "AdaptiveRepulsionData.h"

REGISTER_PARTICLESTUFF(ParticleShaderCropped,"Shader:ParticleShaderCropped");

ParticleShaderCropped::ParticleShaderCropped(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderCropped"))
{

	quad=gluNewQuadric();
	
	new PSParamString(this,&radius_source,"AdaptiveRepulsionData:radius","radattr","radius source",
		"Source of the radius of the form Attribute:param");
	new PSParamDouble(this,&scale,1.35,"scale","scale factor",
		"Amount to scale the radius");
	new PSParamInt(this,&sides,10,"sides","polygon sides",
		"# of sides of the disk approximating polygons");
	new PSParamDouble(this,&toleranceScale,2.0,"tolerance scale",
		"Number of radii we limit fragments to");

	radius_attr = NULL;
	radius_data = NULL;

	pbuffer = NULL;

	// Get the current viewport size
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	// Check to make sure the necessary extensions are supported
	pbuffer->testExtensions();

	// Setup the pbuffer
	setupBuffer(viewport[2], viewport[3]);

	// Setup the depth texture
	glGenTextures(1, &depth_texture);	
	glBindTexture(GL_TEXTURE_RECTANGLE_NV, depth_texture);
	glTexParameteri(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	// Setup the CG environment
	setupCG();

}

void ParticleShaderCropped::attachAttributes()
{
	ParticleShader::attachAttributes();

	std::string attr_name = radius_source.substr(0,radius_source.find(':'));
	std::string param_name = radius_source.substr(radius_source.find(':')+1);
	radius_attr = ps->attributes[attr_name];
	if (radius_attr) {
		int i;
		for (i = 0; i < radius_attr->perparticle.size(); i++) {
			if (radius_attr->perparticle[i]->shortname == param_name) {
				radius_data = (DoubleVector *)radius_attr->perparticle[i]->ref();
				break;
			}
		}
		// If param not found, toss attribute too
		if (i == radius_attr->perparticle.size())
			radius_attr = NULL;
	}

	attachAttribute(imp_int,std::string("ImplicitInterrogator"));
}

/// parameters
int ParticleShaderCropped::qlen()
{
	return 3;
}

void ParticleShaderCropped::getq(double *q)
{
	q[0] = sides;
	q[1] = scale;
	q[2] = toleranceScale;
}

void ParticleShaderCropped::setq(double *q)
{
	sides = (int)q[0];
	scale = q[1];
	toleranceScale = q[2];
}

void ParticleShaderCropped::qname(char **qn)
{
	qn[0] = "Number of sides";
	qn[1] = "Size of disk";
	qn[2] = "Tolerance scale factor";
}

void ParticleShaderCropped::setupBuffer(int width, int height)
{
	hdc = wglGetCurrentDC();
	hglrc = wglGetCurrentContext();

	pbuffer = new nv_pbuffer(width, height);
		
	pbuffer->wglGetLastError();	

}

void ParticleShaderCropped::bindCGParameters()
{

	// Bind the vertex program parameters
	modelViewProj = cgGetNamedParameter(vProgram,"modelViewProj");
	modelView = cgGetNamedParameter(vProgram,"modelView");
	modelViewIT = cgGetNamedParameter(vProgram,"modelViewIT");
	viewMatrix = cgGetNamedParameter(vProgram,"viewMatrix");

	// Bind the fragment program parameters
	lightPositionEC = cgGetNamedParameter(fProgram,"lightPositionEC");
	diffuseColor = cgGetNamedParameter(fProgram,"diffuseColor");
	diskCenterWC = cgGetNamedParameter(vProgram,"diskCenterWC");
	diskRadius = cgGetNamedParameter(fProgram,"diskRadius");
	firstPassDepth = cgGetNamedParameter(fProgram,"firstPassDepth");
	depthRange = cgGetNamedParameter(fProgram,"depthRange");
	toleranceScaleFactor = cgGetNamedParameter(fProgram,"toleranceScaleFactor");

	// Make sure everything got loaded okay
	if (!modelViewProj || !modelView || !lightPositionEC || !diffuseColor || !diskCenterWC 
		|| !diskRadius || !firstPassDepth || !depthRange || !viewMatrix  || ! modelViewIT
		|| !toleranceScaleFactor)
		exit(-1);

	// Setup the texture parameter
	cgGLSetTextureParameter(firstPassDepth, depth_texture);

}

void ParticleShaderCropped::setupCG()
{

	if (cgGLIsProfileSupported(CG_PROFILE_VP30))
		vProfile = CG_PROFILE_VP30;
	else exit(-1);

	if (cgGLIsProfileSupported(CG_PROFILE_FP30))
		fProfile = CG_PROFILE_FP30;
	else exit(-1);

	// Create the Cg context
	CGcontext cgContext = cgCreateContext();
	
	// A VERY common error will occur with these calls: the current directory may easily
	// be changed, causing these two calls to fail.  This will eventually need to be fixed, 
	// but setting the current directory seems to require some windows-specific code, which
	// I am eager to leave out!
	// Create the vertex program
	vProgram = cgCreateProgramFromFile(cgContext, CG_SOURCE, "../Particles/vertex.cg", 
		vProfile, NULL, NULL);
	cgGLLoadProgram(vProgram);

	// Create the fragment program
	fProgram = cgCreateProgramFromFile(cgContext, CG_SOURCE, "../Particles/fragment.cg",
		fProfile, NULL, NULL);
	cgGLLoadProgram(fProgram);

	bindCGParameters();

}

void ParticleShaderCropped::draw(int s)
{

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	// If the window size has changed, we need to reallocate the pbuffer
	if ((pbuffer->width != viewport[2])||(pbuffer->height != viewport[3])) {
		delete pbuffer;
		setupBuffer(viewport[2], viewport[3]);
	}

	// Switch to pbuffer rendering context
	if (wglMakeCurrent(pbuffer->hdc, pbuffer->hglrc) == FALSE)
		pbuffer->wglGetLastError();

	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

	// Draw into the buffer for the depth texture
	for (unsigned int i = 0; i < ps->size(); i++)
		drawParticle(i, s);

	// Switch back to the window rendering context
	if (wglMakeCurrent(hdc, hglrc) == FALSE)
		pbuffer->wglGetLastError();

	// Turn on Cg
	cgGLEnableProfile(vProfile);
	cgGLEnableProfile(fProfile);
	cgGLBindProgram(vProgram);
	cgGLBindProgram(fProgram);

	GLfloat lightpos[4];
	glGetLightfv(GL_LIGHT0, GL_POSITION, lightpos);

	// Set all the Cg parameters that stay constant over all particles
	// Changes to the frustum in DisplayWindow will break this next line
	cgGLSetParameter1f(depthRange, GLfloat(abs(25.0 - 1.0)));
	cgGLSetParameter4fv(lightPositionEC, lightpos);
	cgGLSetParameter1f(toleranceScaleFactor, GLfloat(toleranceScale));
	cgGLSetStateMatrixParameter(viewMatrix, CG_GL_MODELVIEW_MATRIX, CG_GL_MATRIX_IDENTITY);

	// Setup the texture
	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
    glEnable(GL_TEXTURE_RECTANGLE_NV);
	glMatrixMode(GL_MODELVIEW);
	glBindTexture(GL_TEXTURE_RECTANGLE_NV, depth_texture);

	PFNWGLBINDTEXIMAGEARBPROC wglBindTexImageARB =
	  (PFNWGLBINDTEXIMAGEARBPROC)wglGetProcAddress("wglBindTexImageARB");
	
	if (wglBindTexImageARB(pbuffer->hpbuffer, WGL_DEPTH_COMPONENT_NV) == FALSE)
		pbuffer->wglGetLastError();

	cgGLEnableTextureParameter(firstPassDepth);

	// Draw the particles to the screen
	for (unsigned int i = 0; i < ps->size(); i++)
		drawSecondPass(i, s);

	// Disable and unbind the texture 
	cgGLDisableTextureParameter(firstPassDepth);
	glDisable(GL_TEXTURE_RECTANGLE_NV);
	PFNWGLRELEASETEXIMAGEARBPROC wglReleaseTexImageARB =
		(PFNWGLRELEASETEXIMAGEARBPROC)wglGetProcAddress("wglReleaseTexImageARB");
	if (wglReleaseTexImageARB(pbuffer->hpbuffer, WGL_DEPTH_COMPONENT_NV) == false)
		pbuffer->wglGetLastError();

	
	cgGLDisableProfile(vProfile);
    cgGLDisableProfile(fProfile);
}


void ParticleShaderCropped::drawParticle(int i, int s)
{
	// to identify each object 
	glLoadName(i);

	// push the shader name down
	glPushName(0xffffffff);
	glLoadName(s);

	glPushMatrix();

	applyTransformations(i);
	
	drawShape(i);
	
	glPopMatrix();
	glPopName();
}

void ParticleShaderCropped::drawSecondPass(int i, int s)
{

	// to identify each object 
	glLoadName(i);

	// push the shader name down
	glPushName(0xffffffff);
	glLoadName(s);

	glPushMatrix();

	applyTransformations(i);
	setColor(i);

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


void ParticleShaderCropped::drawShape(int i)
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


ParticleShaderCropped::~ParticleShaderCropped()
{
	if (pbuffer) delete pbuffer;
	gluDeleteQuadric(quad);
	cgDestroyContext(cgContext);
}
