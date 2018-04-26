#ifdef WB_USE_CG
#include "ParticleShaderDiskDevelop.h"
#include "AdaptiveRepulsionData.h"

REGISTER_PARTICLESTUFF(ParticleShaderDiskDevelop,"Shader:ParticleShaderDiskDevelop");

ParticleShaderDiskDevelop::ParticleShaderDiskDevelop(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderDiskDevelop"))
{
	quad=gluNewQuadric();

	new PSParamString(this,&radius_source,"AdaptiveRepulsionData:radius","radattr","radius source",
		"Source of the radius of the form Attribute:param");
	new PSParamDouble(this,&scale,0.5,"scale","scale factor",
		"Amount to scale the radius");

	radius_attr = NULL;
	radius_data = NULL;

	setupCG();
}

void ParticleShaderDiskDevelop::attachAttributes()
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

void ParticleShaderDiskDevelop::bindCGParametersVertex()
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


void ParticleShaderDiskDevelop::bindCGParametersFragment()
{

	// Bind the first fragment program parameters
	lightPositionEC = cgGetNamedParameter(fProgram,"lightPositionEC");
	diffuseColor = cgGetNamedParameter(fProgram,"diffuseColor");
	diskRadius = cgGetNamedParameter(fProgram,"diskRadius");
	
	if (!lightPositionEC || !diffuseColor || !diskRadius)
		exit(-1);
}

void ParticleShaderDiskDevelop::setupCG()
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
	fProgram = cgCreateProgramFromFile(cgContext, CG_SOURCE, "../Particles/clip_frag_d.cg",
		fProfile, NULL, NULL);
	cgGLLoadProgram(fProgram);

}

void ParticleShaderDiskDevelop::drawParticle(int i)
{

	cgGLEnableProfile(vProfile);
	cgGLEnableProfile(fProfile);
	cgGLBindProgram(vProgram);
	bindCGParametersVertex();
	bindCGParametersFragment();
	cgGLBindProgram(fProgram);

	GLfloat lightpos[4];
	glGetLightfv(GL_LIGHT0, GL_POSITION, lightpos);

	cgGLSetParameter4fv(lightPositionEC, lightpos);

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
	
	gluDisk(quad,0,radius*scale*1.51,4,4);

	cgGLDisableProfile(vProfile);
    cgGLDisableProfile(fProfile);
}

ParticleShaderDiskDevelop::~ParticleShaderDiskDevelop()
{
	gluDeleteQuadric(quad);
}
#endif