/**
 * Implementation of a shader for silhouette contours
 * @file ShaderSilhouetteContour.cpp
 * @date 6/25/2005
 * @author Matei N. Stroila
 * @remarks
 */
#include "ShaderSilhouetteContour.h"
#include "Attributes/ParticlePosition.h"
#include "Attributes/ParticleVelocity.h"
#include "Attributes/ImplicitInterrogator.h"
#include "Attributes/ViewDependence.h"
#include "Attributes/ParticleVisibility.h"
#include "Attributes/ParticleLocality.h"
#include "Attributes/AdaptiveRepulsionData.h"
#include "Attributes/ParticleVector.h"
#include "Attributes/ParticleBoundingBox.h"



REGISTER_PARTICLESTUFF(ShaderSilhouetteContour,"Shader:ShaderSilhouetteContour");

ShaderSilhouetteContour::ShaderSilhouetteContour(Particles *ps)
:ShaderContour(ps,std::string("ShaderSilhouetteContour"))
{		
	new Attached<ParticleVector>(this,&tangentAttr);
}

void ShaderSilhouetteContour::attachAttributes()
{
	ParticleShader::attachAttributes();
	vis->setOffset(0.5);
	vis->setUseNormalsTest(false);
}

void ShaderSilhouetteContour::findTangent(const unsigned int i)
{
	gmVector3 ni, posi;
	posi = position->getPosition(i);
	ni = impInt->grad(i);
	ni.normalize();
	gmVector3 viewVector(*myCameraPosition -  posi);
	viewVector.normalize();

	gmVector3 tangent = impInt->hess(i) * viewVector; 
	tangent = cross(ni, tangent);

	if(tangent.lengthSquared() > gmEPSILON) tangent.normalize();
	
	tangentAttr->setVector(i,tangent);
}

void ShaderSilhouetteContour::checkCriticalPointPC()
{
	parabolicContoursPar = ps->particleSystem->findParticles(parabContoursName);
	if (!parabolicContoursPar) return;
	ParticlePosition *pcPos = 
parabolicContoursPar->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
	if (!pcPos) return;
	//get the camera position
	gmVector3 cameraPos3 = *(view->getCameraPosition());
	
	gmVector3 gradienti, camerai, xi;
	//SingularityClassification* singClassif;
	for (size_t i = 0; i < parabolicContoursPar->size(); i++)
	{
		xi = pcPos->getPosition(i);
		//ATTENTION these are coming from a different system!!!
		//Therefore we cannot use the precalculated values!
		gradienti = impInt->getImplicit()->grad(xi);
		
		camerai = cameraPos3 - xi;
		double dot_p = fabs(dot(gradienti.normalize(),camerai.normalize()));	  
		if(dot_p< singThreshold ) 
		{
			TNTPoint critPoint(5);
			for(size_t j = 0; j <3; j++){
				critPoint[j] = xi[j];
			}
			critPoint[3] =  0.0; //For now, no Euler angles assigned here
			critPoint[4] =	0.0; //For now, no Euler angles assigned here
			cpList.push_back(ClassifiedPoint(critPoint, UNCLASSIFIED));
			critPointFlag = true;
			//add a silhouette particle at the critical point's position
			this->ps->addParticle();
			position->setPosition(this->ps->size() - 1,xi);
			ps->fullUpdate();
		}
	}

}


void ShaderSilhouetteContour::checkCriticalPoint5D(double *pstatus)
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
	
	p.F = impInt->getImplicit();	;
	p.radius = zoom;

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

	ContourCriticalPointsGSL cpGSL(&p,NewtonTolerance);
	
	numberCP = 0;
	numberCP = cpGSL.FindZeros5IA(theBoxBounds,&p,SILHOUETTES, cpList, pstatus);
	
	if(numberCP > 0)
	{
		printCPs();
		critPointFlag = true;
		ps->fullUpdate();
	}		
	
}

void ShaderSilhouetteContour::checkCriticalPoint4D(double *pstatus)
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
	
	//set parameters
	p.F = impInt->getImplicit();	;
	p.radius = zoom;
		
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
	
	
	ContourCriticalPointsGSL cpGSL(&p,NewtonTolerance);
	numberCP = 0;
	numberCP = cpGSL.FindZeros4D(theBoxBounds,&p,SILHOUETTES, cpList, pstatus);
	if(numberCP > 0)
	{
		printCPs();
		critPointFlag = true;
		ps->fullUpdate();
	}	
}




void ShaderSilhouetteContour::drawPost(void)
{
	
	//draw critical points (if any)
	glDisable(GL_LIGHTING);
	glColor4f((GLfloat) cpColor[0],(GLfloat) cpColor[1],(GLfloat) cpColor[2],(GLfloat) cpColor[3]);
	ClassifiedPointList::const_iterator iter;
	glPointSize(6);
	glBegin(GL_POINTS);
	for(iter = cpList.begin(); iter != cpList.end(); ++iter){	
		glVertex3d((*iter).first[0], (*iter).first[1], (*iter).first[2] );
	}
	glEnd();
	glEnable(GL_LIGHTING);
	if(!drawLines) return;
	//if a critical point is found or the loops connectivity update is required
	if(updateConnectivity || critPointFlag) 
	{
		findChains();
		if(critPointFlag)
		{
			critPointFlag = false;
		}
	}
	
	glLineWidth((GLfloat) lineWidth);
	// draw all lines
	glDisable(GL_LIGHTING);
	glColor4f((GLfloat) color[0],(GLfloat) color[1],(GLfloat) color[2],(GLfloat) color[3]);
	if(!useLoops)
	{
		for(unsigned int i=0;i<(contours->chains).size();i++)
		{
			int isVisble[3] = {1,1,1};
			if(useVisibility) 
			{
				if(cameraPosChanged || !gmIsZero((velocity->v[i]).lengthSquared())) 
				{
					isVisble[0] = vis->getVisibility(myCameraPosition, i);
					isVisble[1] = vis->getVisibility(myCameraPosition, (contours->chains)[i].n1);
					isVisble[2] = vis->getVisibility(myCameraPosition, (contours->chains)[i].n2);
				}				
				else
				{
					isVisble[0] = vis->getVisibility(i);
					isVisble[1] = vis->getVisibility((contours->chains)[i].n1);
					isVisble[2] = vis->getVisibility((contours->chains)[i].n2);
				}
				
			}//end if use visibility
			gmVector3 v1,v2,vi;
			if(isVisble[1] == 1)
				v1 = position->getPosition((contours->chains)[i].n1);
			if(isVisble[2] == 1)
				v2 = position->getPosition((contours->chains)[i].n2);
			if(isVisble[0] == 1)
				vi = position->getPosition(i);
			
			glBegin(GL_LINES);
			
			{
				if(isVisble[0] == 1 && isVisble[1] == 1)
				{
					glVertex3d(v1[0],v1[1],v1[2]);
					glVertex3d(vi[0],vi[1],vi[2]);
				}
				if(isVisble[0] == 1 && isVisble[2] == 1)	
				{
					glVertex3d(v2[0],v2[1],v2[2]);
					glVertex3d(vi[0],vi[1],vi[2]);
				}
			}
			glEnd();
		}//end for 
	}//end if not useLoops
	 //if the connectivity is used draw the loops after checking the visibility 
	 //of the loops
	else
	{
		std::vector<Loop>::const_iterator loops_iter;
		Loop::const_iterator loop_iter, loop_iter_next;
		int isVisble1 =1;
		int isVisble2 =1;
				
		glBegin(GL_LINES);
		for(loops_iter = (contours->loops).begin(); loops_iter != (contours->loops).end(); loops_iter++)
		{
			glColor3f(1,0,0);
			
			loop_iter_next = loops_iter->begin();
			loop_iter_next++;
			for(loop_iter = loops_iter->begin(); loop_iter_next != loops_iter->end(); loop_iter++)
			{
				if(useVisibility) 
				{
					if(cameraPosChanged || !gmIsZero((velocity->v[*loop_iter]).lengthSquared()))
					{
						isVisble1 = vis->getVisibility(myCameraPosition, *loop_iter);
						isVisble2 = vis->getVisibility(myCameraPosition, *loop_iter_next);
					}
					else
					{
						isVisble1 = vis->getVisibility(*loop_iter);
						isVisble2 = vis->getVisibility(*loop_iter_next);
					}
					//end if use visibility	
				}				
				if((isVisble1 == 1) && (isVisble2 == 1))
				{
					gmVector3 v1(position->getPosition(*loop_iter));
					gmVector3 v2(position->getPosition(*loop_iter_next));
					glVertex3d(v1[0],v1[1],v1[2]);
					glVertex3d(v2[0],v2[1],v2[2]);
				}
				loop_iter_next++;
			}//end loop iter

			/*Close loop*/

			if(useVisibility) 
				{
					if(cameraPosChanged)
					{
						isVisble1 = vis->getVisibility(myCameraPosition, (*loops_iter)[loops_iter->size()-1]);
						isVisble2 = vis->getVisibility(myCameraPosition, *(loops_iter->begin()));
					}
					else
					{
						isVisble1 = vis->getVisibility((*loops_iter)[loops_iter->size()-1]);
						isVisble2 = vis->getVisibility(*(loops_iter->begin()));
					}
					//end if use visibility	
				}				
				if((isVisble1 == 1) && (isVisble2 == 1))
				{
					gmVector3 v1(position->getPosition((*loops_iter)[loops_iter->size()-1]));
					gmVector3 v2(position->getPosition( *(loops_iter->begin())));
					glVertex3d(v1[0],v1[1],v1[2]);
					glVertex3d(v2[0],v2[1],v2[2]);
				}

			 //draw arrow to show the loop orientation
			if(!useVisibility && loops_iter->size()>2 )
			{
				gmVector3 v1(position->getPosition(*(loops_iter->begin())));
				gmVector3 v2(position->getPosition((*loops_iter)[2]));
				
				gmVector3 vup = v1 + 0.25 * /*averageLength * */  cross(tangentAttr->getVector(*(loops_iter->begin())),
																  (*myCameraPosition - v1)).normalize();
				gmVector3 vdown = 2 * v1 - vup;
				
				glColor3f(0,1,0);
				glVertex3d(vup[0],vup[1],vup[2]);
				glVertex3d(v2[0],v2[1],v2[2]);
				glVertex3d(vdown[0],vdown[1],vdown[2]);
				glVertex3d(v2[0],v2[1],v2[2]);
			}
		}//end loops iter
		glEnd();
	}//end else if useLoops
	glEnable(GL_LIGHTING);
}

void ShaderSilhouetteContour::checkCriticalPoint_kGkR_4D(double *pstatus)
{
	Box<double> theBoxBounds(4);
	
	//get the implicit bounding box, only once! for a given implicit:
	/*
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
	*/

	//set parameters
	p.F = impInt->getImplicit();
	p.radius = zoom;
	
	//transform Euler angles from degrees to radians
	xrot_last *= gmDEGTORAD;
	yrot_last *= gmDEGTORAD;
	xrot *= gmDEGTORAD;
	yrot *= gmDEGTORAD;
	
	p.phi0 = xrot_last;
	p.phi1 = xrot;
	p.theta0 = yrot_last;
	p.theta1 = yrot;

#if 0	
	//check for already computed CPs
	Intervald Phi(xrot_last, xrot);
	Intervald Theta(yrot_last, yrot);
	ClassifiedPointList::const_iterator iter;
	for(iter = cpList.begin(); iter != cpList.end(); ++iter){	
		if(in((*iter).first[3], Phi) && in((*iter).first[4],Theta))
		{
			std::cout << "found a precomputed CP" << std::endl;
			critPointFlag = true;
			ps->fullUpdate();
			return;
		}
	}
#endif

	theBoxBounds[3].setInterval(0.0,1.0); //time interval
	numberCP = 0;
	ContourCriticalPointsGSL cpGSL(&p,NewtonTolerance);
	for(unsigned int i=0;i<ps->size();i++)
	{
		gmVector3 posi(position->getPosition(i));
		double ri = 3.0 * ardata->r[i];
		for(unsigned int j=0;j<3;j++)
			theBoxBounds[j].setInterval(posi[j]- ri , posi[j] + ri );
		
		numberCP += cpGSL.FindZeros_kGkR_4D(theBoxBounds,&p,SILHOUETTES, cpList, pstatus);
		
	}
	if(numberCP > 0)
	{
		printCPs();
		critPointFlag = true;
		ps->fullUpdate();
	}	
}

void ShaderSilhouetteContour::checkCriticalPointInterval()
{
	Box<double> theBoxBounds(4);
	
	//set parameters
	p.F = impInt->getImplicit();
	p.radius = zoom;
	
	//transform Euler angles from degrees to radians
	xrot_last *= gmDEGTORAD;
	yrot_last *= gmDEGTORAD;
	xrot *= gmDEGTORAD;
	yrot *= gmDEGTORAD;
	
	p.phi0 = xrot_last;
	p.phi1 = xrot;
	p.theta0 = yrot_last;
	p.theta1 = yrot;
#if 0	
	//check for already computed CPs
	Intervald Phi(xrot_last, xrot);
	Intervald Theta(yrot_last, yrot);
	ClassifiedPointList::const_iterator iter;
	for(iter = cpList.begin(); iter != cpList.end(); ++iter){	
		if(in((*iter).first[3], Phi) && in((*iter).first[4],Theta))
		{
			std::cout << "found a precomputed CP" << std::endl;
			critPointFlag = true;
			ps->fullUpdate();
			return;
		}
	}
#endif	
	theBoxBounds[3].setInterval(0.0,1.0); //time interval
	ContourCritical cpCC(&p,NewtonTolerance, useNewton);
	numberCP = 0;
	for(unsigned int i=0;i<ps->size();i++)
	{
		gmVector3 posi(position->getPosition(i));
		double ri = 3.0 * ardata->r[i];
		for(unsigned int j=0;j<3;j++)
			theBoxBounds[j].setInterval(posi[j]- ri , posi[j] + ri );
		
		numberCP += cpCC.search(theBoxBounds,SILHOUETTES, cpList);
	}
	if(numberCP > 0)
	{
		printCPs();
		critPointFlag = true;
		ps->fullUpdate();
	}	
}




