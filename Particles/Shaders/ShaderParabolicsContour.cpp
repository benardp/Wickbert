/**
 * Implementation of a shader for parabolic contours
 * @file ShaderParabolicsContour.cpp
 * @date 12/9/2005
 * @author Matei N. Stroila
 * @remarks
 */

#include "ShaderParabolicsContour.h"
#include "Surface/Implicit/Operator/ZeroGaussianCurvature.h"
#include "Attributes/ParticlePosition.h"
#include "Attributes/ParticleVelocity.h"
#include "Attributes/ImplicitInterrogator.h"
#include "Attributes/ViewDependence.h"
#include "Attributes/ParticleVisibility.h"
#include "Attributes/ParticleLocality.h"
#include "Attributes/AdaptiveRepulsionData.h"
#include "Attributes/ParticleVector.h"
#include "Attributes/ParticleBoundingBox.h"



REGISTER_PARTICLESTUFF(ShaderParabolicsContour,"Shader:ShaderParabolicsContour");

ShaderParabolicsContour::ShaderParabolicsContour(Particles *ps)
:ShaderContour(ps,std::string("ShaderParabolicsContour"))
{		
	imp1 = imp2 = NULL;
	ImplicitIsSet = false;
}

void ShaderParabolicsContour::attachAttributes()
{
	ParticleShader::attachAttributes();
	vis->setOffset(0.5);
	vis->setUseNormalsTest(false);
	ImplicitIsSet = false;
}

void ShaderParabolicsContour::findTangent(const unsigned int i)
{
	if(!ImplicitIsSet) ImplicitIsSet = setImplicit();
	if(!ImplicitIsSet) return;
	gmVector3 n1i, n2i, posi;
	posi = position->getPosition(i);
	n1i = imp1->grad(posi);
	n2i = imp2->grad(posi);	
	gmVector3 tangent;
	tangent = cross(n1i, n2i);
	if(tangent.lengthSquared() > gmEPSILON) tangent.normalize();
	tangentAttr->setVector(i,tangent);
	
}


void ShaderParabolicsContour::drawPost(void)
{
	//draw critical points (if any)
	glDisable(GL_LIGHTING);
	glColor3f(1.0,0.0,0.0);
	ClassifiedPointList::const_iterator iter;
	glPointSize(4);
	glBegin(GL_POINTS);
	for(iter = cpList.begin(); iter != cpList.end(); ++iter){	
		glVertex3d((*iter).first[0], (*iter).first[1], (*iter).first[2] );
	}
	glEnd();
	glEnable(GL_LIGHTING);
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
				
				gmVector3 vup = v1 + 0.25  * cross(tangentAttr->getVector(*(loops_iter->begin())),
																  (*myCameraPosition - v1)).normalize();
				gmVector3 vdown = 2 * v1 - vup;
				
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

bool ShaderParabolicsContour::setImplicit(){
	
	// get the implicit surface from the interrogator
	imp1 = impInt->getImplicit();
	if (!imp1) return false;
	// the second implicit is the corresponding parabolic surface
	if(imp2) delete imp2;
	imp2 = new ZeroGaussianCurvature(imp1);
	return true;
}

ShaderParabolicsContour::~ShaderParabolicsContour()
{
	if(imp2) delete imp2;
}

