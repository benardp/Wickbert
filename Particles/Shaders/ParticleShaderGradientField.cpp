#include "ParticleShader.h"
#include "ImplicitInterrogator.h"
#include "ParticleBoundingBox.h"
#include "ParticleShaderGradientField.h"

REGISTER_PARTICLESTUFF(ParticleShaderGradientField,"Shader:ParticleShaderGradientField");

ParticleShaderGradientField::ParticleShaderGradientField(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderGradientField"))
{
	new PSParamDouble(this,&vectorLength,0.1,"len",
		"Vector length","Length of the gradient glyphs");
	new PSParamDouble(this,&stepSize,1.0,"res","Resolution",
		"Distance between gradient samples.");
	new PSParamBool(this,&normalizeVector,false,"norm","Unitize",
		"Indicate direction but not length of gradients.");

	new Attached<ImplicitInterrogator>(this, &impInt);
	new Attached<ParticleBoundingBox>(this, &bounds);
}

// draw gradient field
void ParticleShaderGradientField::draw(int s)
{
	// check if there is a surface to draw
	if (impInt==NULL || bounds==NULL)
		return;

	gmVector3 &low=bounds->min;
	gmVector3 &high=bounds->max;

	// turn off the lighting
	glPushAttrib(GL_POINT_BIT);
	glDisable(GL_LIGHTING);
	for(double x=low[0];x<=high[0];x+=stepSize)
	for(double y=low[1];y<=high[1];y+=stepSize)
	for(double z=low[2];z<=high[2];z+=stepSize)
	{
		glBegin(GL_LINES);
		gmVector3 origin(x,y,z);
		gmVector3 grad=(impInt->getImplicit()->grad(origin));
		if (normalizeVector)
			grad.normalize();
		grad=grad*vectorLength+origin;
		glVertex3d(origin[0],origin[1],origin[2]);
		glVertex3d(grad[0],grad[1],grad[2]);
		glEnd();
		// a point as the orgin of
		glPointSize(4);
		glBegin(GL_POINTS);
		glVertex3d(origin[0],origin[1],origin[2]);
		glEnd();
	}
	glEnable(GL_LIGHTING);
	glPopAttrib();
}