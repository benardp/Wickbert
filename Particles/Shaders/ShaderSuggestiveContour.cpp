/*
 *  ShaderSuggestiveContour.cpp
 *  gps
 *
 *  Created by Matei Stroila on 7/16/05.
 * 
 *
 */
#include <vector>
#include <set>
#include "ShaderSuggestiveContour.h"
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


REGISTER_PARTICLESTUFF(ShaderSuggestiveContour,"Shader:ShaderSuggestiveContour");

ShaderSuggestiveContour::ShaderSuggestiveContour(Particles *ps)
:ShaderContour(ps,std::string("ShaderSuggestiveContour"))
{
	new Attached<LightPosition>(this,&light);
	new Attached<ParticleVelocity>(this,&velocity);
	new PSParamDouble(this,&velocityThreshold,100.0,"velocityThreshold","Velocity Threshold");
	new PSParamDouble(this,&dirDerivThreshold,-0.01,"Threshold","Dir Deriv Threshold");
}


void ShaderSuggestiveContour::drawPost(void)
{
	//if unstable, do not draw
	for(unsigned int i = 0; i < ps->size(); i++)
	{
		if((velocity->v[i]).lengthSquared() > velocityThreshold) return;
	}
	
	double dirDeriv1, dirDeriv2, dirDerivi; //directional derivatives
	if(updateConnectivity || critPointFlag) 
	{
		findChains();
		if(critPointFlag)
		{
			critPointFlag = false;
			return;
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
			dirDeriv1 = dirDeriv2 = dirDerivi = 1.0;
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
				if(isVisble[1] == 1)
				{
					dirDeriv1 = computeDerivative((contours->chains)[i].n1);
				}
				if(isVisble[2] == 1)
				{
					dirDeriv2 = computeDerivative((contours->chains)[i].n2);
				}
				if(isVisble[0] == 1)
				{
					dirDerivi = computeDerivative(i);
				}
			}//end if use visibility
			gmVector3 v1,v2,vi;
			
			if(isVisble[1] == 1)
			{
				v1 = position->getPosition((contours->chains)[i].n1);
			}
			if(isVisble[2] == 1)
			{
				v2 = position->getPosition((contours->chains)[i].n2);
			}
			if(isVisble[0] == 1)
			{
				vi = position->getPosition(i);
			}
			
			glBegin(GL_LINES);
			
			{
				if(isVisble[0] == 1 && isVisble[1] == 1 && 
				   dirDerivi>0 && dirDeriv1>0)
				{
					glVertex3d(v1[0],v1[1],v1[2]);
					glVertex3d(vi[0],vi[1],vi[2]);
				}
				if(isVisble[0] == 1 && isVisble[2] == 1 && 
				   dirDerivi>0 && dirDeriv2>0)	
				{
					glVertex3d(v2[0],v2[1],v2[2]);
					glVertex3d(vi[0],vi[1],vi[2]);
				}
			}
			glEnd();
		}//end for 
	}//end if not useLoops
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
					if(isVisble1 == 1)
					{
						if(computeDerivative(*loop_iter) < dirDerivThreshold)
						{
							vis->setVisibility(*loop_iter,0);
							isVisble1 = 0;
						}
					}
					if(isVisble2 == 1)
					{
						if(computeDerivative(*loop_iter_next) < dirDerivThreshold)
						{
							vis->setVisibility(*loop_iter_next,0);
							isVisble2 = 0;
						}
					}
				}	//end if useVisibility
				
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

			//draw arrows
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


void ShaderSuggestiveContour::findTangent(const unsigned int i)
{
	gmVector3 ni, posi;
	posi = position->getPosition(i);
	ni = impInt->grad(i);
	ni.normalize();

	gmVector3 v = *(view->getCameraPosition()) - posi; //view vector
	gmMatrix3 Hess = impInt->hess(i); //hessian
	gmVector3 Hessv = Hess * v;
	Implicit* imp = impInt->getImplicit();
	gmVector3 ugly(
				   dot((imp->hessi(posi,0) * v),v),
				   dot((imp->hessi(posi,1) * v),v),
				   dot((imp->hessi(posi,2) * v),v)
				   );

	gmVector3 tangent =  ugly - 2 * Hessv ;
	tangent.normalize();
	tangent = cross(ni,tangent);

	if(tangent.lengthSquared() > gmEPSILON) tangent.normalize();	
	tangentAttr->setVector(i,tangent);

}


double ShaderSuggestiveContour::computeDerivative(const int i)
{
	double dirDeriv;	//directional derivative (to be computed)
	gmVector3 v;	//view direction
	gmVector3 x = position->getPosition(i);	//particle position
	gmVector3 suggNormal; //gradient of the suggestive surface
	
	v = *myCameraPosition - x; //view vector
	
	gmMatrix3 Hess = impInt->hess(i);//getImplicit()->hess(x) ; //hessian
	gmVector3 Hessv = Hess * v;
	gmVector3 ugly(
				   dot((impInt->getImplicit()->hessi(x,0) * v),v),
				   dot((impInt->getImplicit()->hessi(x,1) * v),v),
				   dot((impInt->getImplicit()->hessi(x,2) * v),v)
				   );
	
	suggNormal = ugly - 2 * Hessv ;
	dirDeriv = dot(suggNormal,v);
	return dirDeriv;
}
