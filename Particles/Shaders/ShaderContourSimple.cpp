/**
 * Implementation of a contour shader 
 * @file ShaderContour.cpp
 * @date 09/08/2005
 * @author Matei N. Stroila
 * @remarks
 */
#include "ShaderContourSimple.h"

#include "Attributes/ParticlePosition.h"
#include "Attributes/ParticleVelocity.h"
#include "Attributes/ImplicitInterrogator.h"
#include "Attributes/ViewDependence.h"
#include "Attributes/ParticleVisibility.h"
#include "Attributes/ParticleLocality.h"
#include "Attributes/AdaptiveRepulsionData.h"
#include "Attributes/ParticleVector.h"
#include "Attributes/ParticleBoundingBox.h"


#include "ism/ContourCriticalPointsGSL.h"
#include "ism/ContourCritical.h"
	
#include <vector>
#include <set>
#include <queue>



ShaderContourSimple::LoopCallback::LoopCallback (ShaderContourSimple::LoopMethod * lm)	
:PSParamComboBox::Callback()
{
	_lm=lm;
	*_lm=ShaderContourSimple::SIMPLE_LOOPS;
	_choices.push_back("SimpleLoops");
	_choices.push_back("HeapLoops");
}

void ShaderContourSimple::LoopCallback::itemselected(unsigned int selection)
{	
	if (selection==0) //simpleLoops
		*_lm=ShaderContourSimple::SIMPLE_LOOPS;
	else 
		*_lm=ShaderContourSimple::HEAP_LOOPS;
}


ShaderContourSimple::ShaderContourSimple(Particles *ps, const std::string& name)
:ParticleShader(ps,name)
{
	new PSParamgmVector4(this,&color,gmVector4(1.0,0.0,0.0,1.0),
						 "color","Color","Line color");
	new PSParamgmVector4(this,&invColor,gmVector4(1.0,0.0,0.0,1.0),
						 "invColor","Invisible Pen Color","Invisible Line Color");
	new PSParamgmVector3(this,&brushColorf,gmVector3(1.0,0.8,0.6),
						 "brushColorf","Brush color","Brush color");
	new PSParamInt(this,&penWidth,4,"penWidth","Pen width","Pen width");
	new PSParamInt(this,&invPenWidth,0,"invPenWidth","Invisible Pen width","Invisible Pen width");
	new PSParamBool(this,&isShaded,true,"isShaded","Shade","Set true to shade it.");
	
	new PSParamBool(this,&drawLines,true,"drawLines","Draw Lines","Set true to draw GL lines.");
	new PSParamBool(this,&useVisibility,false,"Visibility","Use Visibility","Set true for Visibility check.");
	new PSParamBool(this,&useLoops,true,"useLoops","Use Loops","Set true to build Loops.");
	new PSParamBool(this,&updateConnectivity,true,"updateConnectivity","Update Connectivity","Set true if want connectivity update");
	new PSParamBool(this,&useLoopsCGAL,false,"useLoopsCGAL","Use CGAL Loops","Set true to use CGAL Loops visibility.");

	new PSParamComboBox(this,new LoopCallback(&loopMethod),"loopCall","Loop:","Method to select the way loops are chained.");

	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
	new Attached<ImplicitInterrogator>(this,&impInt);
	new Attached<ParticleVisibility>(this,&vis);
	new Attached<ParticleLocality>(this,&p_locality);
	new Attached<AdaptiveRepulsionData>(this,&ardata);
	new Attached<ViewDependence>(this,&view);
	new Attached<ParticleVector>(this,&tangentAttr);
	new Attached<ParticleBoundingBox>(this,&pBoundingBox);
	new Attached<Contours>(this,&contours);

	lineWidth=2;
	cameraPosChanged = true;
	oldCameraPosition = new gmVector3(0.0,0.0,0.0);
}

void ShaderContourSimple::drawPre()
{
	
	//check whether the camera position has changed
	//if it has changed, recompute the tangent vectors to the silhouette at 
	//all particle postitions
	
	myCameraPosition = view->getCameraPosition();
	*oldCameraPosition = *myCameraPosition;
	cameraPosChanged = true;
	for(unsigned int i=0;i<ps->size();i++)
	{
		findTangent(i);
	}
	return;
}

/// Draws contours. Implemented here in drawPost because it isn't a per-particle
/// operation.
void ShaderContourSimple::drawPost(void)
{
	if(updateConnectivity) 
		findChains();
	
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
			glColor4f(1.0f - (GLfloat) color[0],1.0f - (GLfloat) color[1],1.0f - (GLfloat) color[2],(GLfloat) color[3]);
	
			if(!useVisibility && loops_iter->size()>2 )
			{
				gmVector3 v1(position->getPosition(*(loops_iter->begin())));
				gmVector3 v2(position->getPosition((*loops_iter)[2]));
				
				gmVector3 vup = v1 + 0.25 * cross(tangentAttr->getVector(*(loops_iter->begin())),
																  (*myCameraPosition - v1)).normalize();
				gmVector3 vdown = 2 * v1 - vup;
				
				glVertex3d(vup[0],vup[1],vup[2]);
				glVertex3d(v2[0],v2[1],v2[2]);
				glVertex3d(vdown[0],vdown[1],vdown[2]);
				glVertex3d(v2[0],v2[1],v2[2]);
			}
			glColor4f((GLfloat) color[0],(GLfloat) color[1],(GLfloat) color[2],(GLfloat) color[3]);
		}//end loops iter
		glEnd();
	}//end else if useLoops
	glEnable(GL_LIGHTING);
}

namespace{
class ParticleConnection{
public:	
	ParticleConnection(unsigned int i , unsigned int j, double d)
		:con(i,j)
		,distance(d)
	{};

	std::pair<unsigned int, unsigned int> con;
	double distance;

	bool operator< (const ParticleConnection& right) const 
	{
		return (distance>right.distance); 
	}
};
}
void ShaderContourSimple::findChains()
{
	unsigned int i,j,k;
		
	if (loopMethod==HEAP_LOOPS)
	{
		//we want a special distance function based on the normal and gradient at the point.
		//for the moment we will simplify things and only test for the normal of the implicit surfaces.
		//For specs this is not sufficient, but for the others, this should already do a good job...


		const unsigned int psSize = ps->size();
		std::vector<double> squareDistances;
		//yes, it is potentially too big
		//squareDistances.resize(psSize *(psSize +1)/2,DBL_MAX);
	
		std::priority_queue<ParticleConnection> possibleConnections;
		//we precacalculate all distances
		for(i=0;i<ps->size();i++)
		{
				//if (fabs(velocity->v[i].lengthSquared())>pow(0.01*this->ardata->sigma_max,2))
				//	continue;

				//all neighbors of i in the queryRadius are in neighbors
				std::list<unsigned int> neighbors;
				double queryRadius = ardata->sdmul * ardata->r[i];
				//queryRadius= queryRadius>ardata->sigma_max ? ardata->sigma_max : queryRadius;
				//DEBUG ELMAR we should add a function to the locality which delivers neighbors with squaredDistance,
				//otherwise we have to calculate this twice... because distance is calculated in the locality anyway!
				
				p_locality->getNeighbors(i,queryRadius,neighbors);
				const std::list<unsigned int>::iterator end=neighbors.end();
				for (std::list<unsigned int>::iterator nbr=neighbors.begin();nbr!=end;++nbr)
				{
					//we want all distances only once!
					if (i<*nbr)
					{
						//if (fabs(velocity->v[*nbr].lengthSquared())>pow(0.01*this->ardata->sigma_max,2))
						//	continue;
						
						//now we want to calculate a better distance function.
						//basically we would like to push according to the distance from both surfaces...
						
						const gmVector3 normal = impInt->normal(i);
						gmVector3 diff = position->getPosition(*nbr)-position->getPosition(i);
						possibleConnections.push(ParticleConnection(i,*nbr, diff.lengthSquared()));
						
					}
				}
		}

		//now we chain the closest particles
		std::vector<std::vector<unsigned int> > neighbors;
		neighbors.resize(psSize);

		for (;!possibleConnections.empty();possibleConnections.pop())
		{
			const ParticleConnection & currCon=possibleConnections.top();
			unsigned int firstParticle=currCon.con.first;
			unsigned int secondParticle=currCon.con.second;
	

			if ((neighbors[firstParticle].size()<=2)&&(neighbors[secondParticle].size()<=2))
			{
				neighbors[firstParticle].push_back(secondParticle);
				neighbors[secondParticle].push_back(firstParticle);
			}
		}

		for (i=0;i<psSize;++i)
		{
			Neighbors n;
			if (neighbors[i].size()>1)
			{
				n.n1=(neighbors[i])[0];
				n.n2=(neighbors[i])[1];


				//orient and define chains
				gmVector3 v1,v2;
				v1 = position->getPosition(n.n1);
				v2 = position->getPosition(n.n2);
					
				if((dot(v2 - v1, tangentAttr->getVector(i)) > 0)) //good orientation or no neighbors
					contours->setChains(i, n);
				else
				{
					unsigned int temp = n.n1; //swap n1 and n2 to get good orientation
					n.n1 = n.n2;
					n.n2 = temp;
					contours->setChains(i, n);
				}	
			}
			else
			{
				n.n1=i;
				n.n2=i;

				contours->setChains(i,n);
			}
		}
	}
	else//if (loopMethod==SIMPLE_LOOPS)
	{
		for(i=0;i<ps->size();i++)
		{
			//the two closest neighbors of i are hold by n; need to find them
			Neighbors n;
			n.n1=i;
			n.n2=i;
			
			//all neighbors of i in the queryRadius are in neighbors
			std::list<unsigned int> neighbors;
			std::list<unsigned int>::iterator nbr;
			
			neighbors.clear();

			double queryRadius = ardata->sdmul * ardata->r[i];
			queryRadius= queryRadius>ardata->sigma_max ? ardata->sigma_max : queryRadius;
			p_locality->getNeighbors(i,queryRadius,neighbors);
		
			//much faster than sorting!
			if(neighbors.size()>1)
			{
				unsigned int n1=0;
				unsigned int n2=n1;
				double minN1=DBL_MAX;
				double minN2=DBL_MAX;

				for (nbr = neighbors.begin(); nbr != neighbors.end(); nbr++) 
				{
					j = *nbr;
					gmVector3 v1,v2, delta;
					v1 = position->getPosition(i);
					v2 = position->getPosition(j);
					double tangentWeight = dot((v2 - v1).normalize(), tangentAttr->getVector(i));
					double currDistance=distanceSquared(position->getPosition(j),position->getPosition(i));
					if( tangentWeight > 0){ 
						if (currDistance<minN2){
							n2 = j;
							minN2=currDistance;
						}
					} 
					else { 
						if (currDistance<minN1){
							n1 = j;
							minN1=currDistance;
						}
					} 
				}

				n.n1=n1;
				n.n2=n2;
			}

			contours->setChains(i, n);
		}
	}


	if(useLoops)
	{	
		(contours->loops).clear();
		std::set<int> visited_chains; //memorize the visited particles
		i = 0;
		while(i<ps->size())
		{
			//if already visited increase i else start a loop at i
			if(visited_chains.find(i) != visited_chains.end()) i++;
			else
			{
				//if i is not connected to anything  mark i as visited (i.e. insert it in
				//visited_chains, and go to next i
				//else start the loop
				if(((contours->chains)[i].n1 == i) && ((contours->chains)[i].n2 == i) )
				{
					visited_chains.insert(i);
					i++;
				}
				else
				{
					Loop loop;
					loop.push_back(i);
					visited_chains.insert(i);
					k = (contours->chains)[i].n2; //move to next particle
					loop.push_back(k);
					while((visited_chains.find(k) == visited_chains.end()))
					{
						visited_chains.insert(k);
						k = (contours->chains)[k].n2;
						loop.push_back(k);
					}//end while
					loop.pop_back();
					//correct chain information
					//the neighbor n1 of the first entry will be the last in the loop
					contours->chains[i].n1=loop[loop.size()-1];
					contours->chains[loop[loop.size()-1]].n2=i;
					(contours->loops).push_back(loop);
				}//end if
			}//end if
		}//end while
	}//end if useLoops
}

void ShaderContourSimple::event(int e) 
{ 
}

ShaderContourSimple::~ShaderContourSimple()
{
	delete oldCameraPosition;
}

/** Called after particles are added to the system.
	*/ 
void ShaderContourSimple::particleAdded(){

	if(!updateConnectivity) findChains();
}

	/** Called after particle are removed from the system.
	*/
void ShaderContourSimple::particleRemoved(unsigned int i){

	if(!updateConnectivity) findChains();
}

