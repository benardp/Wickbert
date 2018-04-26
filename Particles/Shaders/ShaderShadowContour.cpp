/**
 * Implementation of a shadow contour shader 
 * @file ShaderShadowContour.h
 * @date 06/25/2005
 * @author Matei N. Stroila
 * @remarks
 */

#include "ShaderShadowContour.h"
#include "Attributes/ParticlePosition.h"
#include "Attributes/ParticleVelocity.h"
#include "Attributes/ImplicitInterrogator.h"
#include "Attributes/ViewDependence.h"
#include "Attributes/ParticleVisibility.h"
#include "Attributes/ParticleLocality.h"
#include "Attributes/AdaptiveRepulsionData.h"
#include "Attributes/ParticleVector.h"
#include "Attributes/ParticleBoundingBox.h"
#include "Attributes/LightPosition.h"



REGISTER_PARTICLESTUFF(ShaderShadowContour,"Shader:ShaderShadowContour");

ShaderShadowContour::ShaderShadowContour(Particles *ps)
:ShaderContour(ps,std::string("ShaderShadowContour"))
{
	new PSParamgmVector3(this,&shadeColorf,gmVector3(0.8,0.6,0.4),
						 "shadeColorf","Shade color","Shade color");
	new Attached<LightPosition>(this,&light);
}

void ShaderShadowContour::findTangent(const unsigned int i)
{
	gmVector3 ni, posi;
	posi = position->getPosition(i);
	ni = impInt->grad(i);
	
	gmVector4 lightPos4;
	light->getLightPosition(lightPos4);
	gmVector3 lightPos3(lightPos4[0],lightPos4[1],lightPos4[2]);	
	
	gmVector3 lightVector;
	gmVector3 tangent;
	if(lightPos4[3]==0)
	{
		//directional
		lightVector = lightPos3;

		tangent = impInt->hess(i) * lightPos3;
		lightVector.normalize();
	}
	else
	{
		lightVector = lightPos3 - lightPos4[3]*posi;
		tangent = impInt->hess(i) * lightVector;
		lightVector.normalize();
	}
	
	tangent.normalize();
	ni.normalize();
	
	tangent = (-1) * cross(ni,tangent);
	
	if(tangent.lengthSquared() > gmEPSILON) tangent.normalize();
	tangentAttr->setVector(i,tangent);
}

void ShaderShadowContour::checkCriticalPointPC()
{
	parabolicContoursPar = ps->particleSystem->findParticles(parabContoursName);
	if (!parabolicContoursPar) return;
	ParticlePosition *pcPos =
		parabolicContoursPar->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
	if (!pcPos) return;
	//get the light position in eye coordinates
	gmVector4 lightPos4;
	light->getLightPosition(lightPos4);
	gmVector3 lightPos3(lightPos4[0],lightPos4[1],lightPos4[2]);	

	
	gmVector3 gradienti, lighti, xi;
	//SingularityClassification* singClassif;
	for (size_t i = 0; i < parabolicContoursPar->size(); i++)
	{
		xi = pcPos->getPosition(i);
		//these are coming from a different system...
		//therefore we cannot use the cache - Elmar
		gradienti = impInt->getImplicit()->grad(xi);
		if(lightPos4[3]==0)
			lighti = lightPos3;
		else
			lighti = lightPos3 - lightPos4[3]*xi;
			  
		if(fabs(dot(gradienti.normalize(),lighti.normalize())) < singThreshold ) 
		{
			TNTPoint critPoint(5);
			for(size_t j = 0; j <3; j++){
				critPoint[j] = xi[j];
			}
			critPoint[3] =  0.0; //For now, no Euler angles assigned here
			critPoint[4] =	0.0; //For now, no Euler angles assigned here
			cpList.push_back(ClassifiedPoint(critPoint, UNCLASSIFIED));
			critPointFlag = true;
			ps->fullUpdate();
		}
	}

}

void ShaderShadowContour::checkCriticalPoint4D(double *pstatus)
{
	static Box<double> theBoxBounds(4);
	
	//get the implicit bounding box, only once! for a given implicit:
	static bool init = true;
	if(init)
	{
		gmVector3 _min,_max;
		pBoundingBox->computeBoundingBox(_min,_max);
		
		for (unsigned int i = 0; i < 3; i++)
			//enlarge the bounding box with a tenth of the bounds
			theBoxBounds[i].setInterval(_min[i] - 0.1 * fabs(_min[i]) ,_max[i] + 0.1 * fabs(_max[i]));
		init = false;
	}
	
	//get the light position in eye coordinates
	GLfloat lightPositionf[4];
	glGetLightfv(GL_LIGHT0, GL_POSITION, lightPositionf);
	
	if(!gmIsZero(lightPositionf[3]))  
	{	
		std::cout << "not a directional light! The positional light is not yet implemented" << std::endl;
		return;
	}
	p.F = impInt->getImplicit();	;
	p.radius = zoom;
	p.Lx = lightPositionf[0];
	p.Ly = lightPositionf[1];
	p.Lz = lightPositionf[2];
	
	//transform Euler angles from degrees to radians
	xrot_last *= gmDEGTORAD;
	yrot_last *= gmDEGTORAD;
	xrot *= gmDEGTORAD;
	yrot *= gmDEGTORAD;
	
	p.phi0 = xrot_last;
	p.phi1 = xrot;
	p.theta0 = yrot_last;
	p.theta1 = yrot;
	theBoxBounds[3].setInterval(0.0,1.0);
	
	ContourCriticalPointsGSL cpGSL(&p, NewtonTolerance);

	numberCP = 0;
	numberCP = cpGSL.FindZeros4D(theBoxBounds,&p,SHADOWS, cpList, pstatus);
	
	if(numberCP > 0)
	{
		printCPs();
		critPointFlag = true;
		ps->fullUpdate();
	}	
	
}

void ShaderShadowContour::checkCriticalPoint5D(double *pstatus)
{
	static Box<double> theBoxBounds(5);
	
	//get the implicit bounding box, only once! for a given implicit:
	static bool init = true;
	if(init)
	{
		gmVector3 _min,_max;
		pBoundingBox->computeBoundingBox(_min,_max);
		
		for (unsigned int i = 0; i < 3; i++)
			//enlarge the bounding box with a tenth of the bounds
			theBoxBounds[i].setInterval(_min[i] - 0.1 * fabs(_min[i]) ,_max[i] + 0.1 * fabs(_max[i]));
		init = false;
	}
	
	//get the light position in eye coordinates
	GLfloat lightPositionf[4];
	glGetLightfv(GL_LIGHT0, GL_POSITION, lightPositionf);
	
	if(!gmIsZero(lightPositionf[3]))  
	{	
		std::cout << "not a directional light! The positional light is not yet implemented" << std::endl;
		return;
	}
	p.F = impInt->getImplicit();	;
	p.radius = zoom;
	p.Lx = lightPositionf[0];
	p.Ly = lightPositionf[1];
	p.Lz = lightPositionf[2];
	
	//transform Euler angles from degrees to radians
	xrot_last *= gmDEGTORAD;
	yrot_last *= gmDEGTORAD;
	xrot *= gmDEGTORAD;
	yrot *= gmDEGTORAD;
	
	//check interval bounds
	if(xrot_last <= xrot)	
		theBoxBounds[3].setInterval(xrot_last,xrot);
	else
		theBoxBounds[3].setInterval(xrot,xrot_last);
	
	if(yrot_last <= yrot)	
		theBoxBounds[4].setInterval(yrot_last,yrot);
	else
		theBoxBounds[4].setInterval(yrot,yrot_last);
	
	std::cout << "the box bounds for the Newton search:" << std::endl;
	std::cout << xrot_last << "," << xrot << "," << yrot_last << "," <<  yrot << std::endl;
	
	ContourCriticalPointsGSL cpGSL(&p, NewtonTolerance);

	numberCP = 0;
	numberCP = cpGSL.FindZeros5D(theBoxBounds,&p,SHADOWS, cpList, pstatus);
	
	if(numberCP > 0)
	{
		printCPs();
		critPointFlag = true;
		ps->fullUpdate();
	}		
	
}



ShaderShadowContour::~ShaderShadowContour(){}
