#include "FindFirstIntersectionWithImplicitsShader.h"
#include "ParticleVisibility.h"
#include "ViewDependence.h"
#include "ImplicitInterrogator.h"
#include "AdaptiveRepulsionData.h"
#include "ParticlePosition.h"

#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>



REGISTER_PARTICLESTUFF(FindFirstIntersectionWithImplicits,"Shader:FindFirstIntersectionWithImplicits");




double FindFirstIntersectionWithImplicits::ImplicitOnLine::proc(double t, void * params)
{
	return imp->proc(start+diff*t);
}

gmVector3 FindFirstIntersectionWithImplicits::ImplicitOnLine::start=gmVector3(0,0,0);
gmVector3 FindFirstIntersectionWithImplicits::ImplicitOnLine::diff=gmVector3(0,0,0);
Implicit * FindFirstIntersectionWithImplicits::ImplicitOnLine::imp;






FindFirstIntersectionWithImplicits::FindFirstIntersectionWithImplicits(Particles *ps)
	:ParticleShader(ps,std::string("FindFirstIntersectionWithImplicits"))
{
	if (ps!=0)
		std::cout<<"Added the find intersection shader to "<<ps->name<<'\n';
	new Attached<ViewDependence>(this,&_view);
	new Attached<ImplicitInterrogator>(this,&_impInt);
	new Attached<ParticlePosition>(this,&_position);

	new PSParamDouble(this,&_stepsPerUnit, 10.0,gmEPSILON,1000000000.0,"step","Step per unit:","steps per unit distance");
	new PSParamInt(this,&_recursive, 5,"recursive","recursive refine:",0,100000000,"recursive refinement in step");
	new PSParamInt(this,&_stepsPerPrecisionTest, 5,"#StepsPerPrecisionTest","#StepsPerPrecisionTest:",0,100000000,"Evaluating for precision is costly, therefore" 
																											"this value specifies how many iterations are made"
																											"before the current value is checked against the"
																											"wanted precision (0 means never)");
	new PSParamDouble(this,&_precision,gmEPSILON,"precision","Precision:","Wanted precision to solve for, it still respects the maximum number of steps.");

	new PSParamgmVector3(this,&_destination,gmVector3(0,0,0),"dest","dest:","destination of ray");
	new PSParamBool(this,&_shootTestRay,false,"testRay","TestRay","shoots a test ray from the eye to the destination (fix viewpoint to turn around and see it)");
	new PSParamBool(this,&_shootTestRaysThroughParticles,false,"testRays","ParticleRays","shoot a test ray through each particle");
}



void FindFirstIntersectionWithImplicits::attachAttributes()
{
	ParticleShader::attachAttributes();
}

void FindFirstIntersectionWithImplicits::drawPre()
{
	gmVector3 * cameraPosition = _view->getCameraPosition();

	if (_shootTestRay)
	{
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glPointSize(5);
		glDisable(GL_LIGHTING);
		gmVector3 diff=_destination-(*cameraPosition);
		double length=diff.length();
		length*=_stepsPerUnit;
		unsigned int steps=(unsigned int) fabs(length);

		Implicit * imp=_impInt->getImplicit();
		std::vector<Implicit*> imps;
		imps.push_back(imp);

		std::vector<Implicit *> intersected;
		gmVector3 positive;
		gmVector3 negative;

		findNegativeStepValue(*cameraPosition,_destination,imps,_stepsPerUnit,&positive, &negative,&intersected);

		if (!intersected.empty())
		{
			std::pair< gmVector3 , Implicit* > result;
			findIntersection(*cameraPosition, _destination,imps,&result, _stepsPerUnit, _recursive, _precision, (_stepsPerPrecisionTest>0)?_stepsPerPrecisionTest:0);

			glColor3f(0.0f,0.0f,0.0f);
			glPointSize(10);
			glBegin(GL_POINTS);
			glVertex3d(positive[0],positive[1],positive[2]);
			glVertex3d(negative[0],negative[1],negative[2]);
			glVertex3d(result.first[0], result.first[1], result.first[2]);
			glEnd();


			glBegin(GL_LINES);
			glColor3f(0.5,1,0.5);
			glVertex3d((*cameraPosition)[0], (*cameraPosition)[1], (*cameraPosition)[2]);
			glVertex3d(result.first[0],result.first[1],result.first[2]);
			
			glColor3f(1,0.5,0.5);
			glVertex3d(result.first[0],result.first[1],result.first[2]);		
			glVertex3d(_destination[0], _destination[1], _destination[2]);		
			glEnd();

		}
		glPopAttrib();
	}
}


bool FindFirstIntersectionWithImplicits::findIntersection(	const gmVector3 & start, 
															const gmVector3 & stop, 
															Implicit* implicit,
															gmVector3 * result,
															double stepsPerUnit, 
															unsigned int recursive, 
															double precision,
															unsigned int stepsPerPrecisionTest)
{
	//a little waste of time, but really little... but more secure and less code copy paste
	assert(implicit);
	std::vector<Implicit*> implicits;
	implicits.push_back(implicit);
	std::pair<gmVector3, Implicit*> res;
	findIntersection(start, stop, implicits,&res,stepsPerUnit,recursive,precision,stepsPerPrecisionTest);

	if (res.second)
	{
		*result=res.first;
		return true;
	}
	else
	{
		return false;
	}
}





void FindFirstIntersectionWithImplicits::findIntersection(	const gmVector3 & start, 
															const gmVector3 & stop, 
															std::vector<Implicit*> implicits,
															std::pair<gmVector3, Implicit*> * result,
															double stepsPerUnit, 
															unsigned int recursive, 
															double precision,
															unsigned int stepsPerPrecisionTest)
{
	assert(result);
	assert(((stepsPerPrecisionTest>0)||(precision>0))||(stepsPerPrecisionTest=0));
	assert(stepsPerUnit>0);
	
	result->second=0;

	gmVector3 positive;
	gmVector3 negative;
	std::vector<Implicit *> intersectedImplicits;
	//we find the first position where the step changes from positive to negative
	findNegativeStepValue(start, stop, implicits, stepsPerUnit, &positive, &negative, &intersectedImplicits);

	//if there is none, we did not hit the surface and we return a 0 pointer
	if (intersectedImplicits.empty())
	{
		result->second=0;
		return;
	}
	//now we found all potential intersection
	//we will now recursively adapt to find the real first intersection

	double minDist=1.0;
	for (std::vector<Implicit *>::iterator iter=intersectedImplicits.begin(); iter!=intersectedImplicits.end();++iter)
	{
		double dist=recursiveRootSearch(positive, negative, *iter, recursive, precision, stepsPerPrecisionTest);
		if (dist<minDist)
		{
			result->second=*iter;
			minDist=dist;
		}
	}
	result->first=positive+(negative-positive)*minDist;
}

void FindFirstIntersectionWithImplicits::findNegativeStepValue(	const gmVector3 & start, 
																const gmVector3 & stop, 
																const std::vector<Implicit *> implicits, 
																double stepsPerUnit, 
																gmVector3 * positionPositive, gmVector3 * positionNegative, 
																std::vector<Implicit*> * intersectedImplicits)
{
	assert(positionPositive);
	assert(positionNegative);
	assert(intersectedImplicits);
	assert(stepsPerUnit>0);

	gmVector3 diff= stop - start;
	double length = diff.length();
	length *= stepsPerUnit;
	unsigned int steps=(unsigned int) fabs(length);	

	diff/=(double)stepsPerUnit;
	*positionNegative=start;
	*positionPositive=start;
	
	
	for (unsigned int i=0;i<steps+1;++i)
	{
		//incremental might lead to less precision...
		*positionNegative+=diff;
		for (std::vector<Implicit *>::const_iterator iter=implicits.begin(); iter!=implicits.end();++iter)
		{
			if ((*iter)->proc(*positionNegative)<0)
			{
				intersectedImplicits->push_back(*iter);
			}
		}
		if (!(intersectedImplicits->empty()))
			break;

		*positionPositive=*positionNegative;
	}
}

double FindFirstIntersectionWithImplicits::recursiveRootSearch(const gmVector3 & start, const gmVector3 & stop , 
															   Implicit * implicit, 
															   unsigned int nbSteps, double precision, unsigned int stepsPerPrecisionTest)
{
	assert(implicit);

	{
	ImplicitOnLine::start=start;
	ImplicitOnLine::diff= stop - start;
	ImplicitOnLine::imp=implicit;


	double r=0;
	double x_hi=1;
	double x_lo=0;
	unsigned int count=0;
	unsigned int iter=0;

	do
	{
		//this seems safer
		if (x_hi-x_lo<gmEPSILON)
			break;
		++iter;
		++count;
		r=(x_hi+x_lo)/2.0;
		double v=ImplicitOnLine::proc(r,0);
		if (v<-gmEPSILON)
			x_hi=r;
		else if (v>gmEPSILON)
			x_lo=r;
		else
			break;
			
		if (count==stepsPerPrecisionTest)
		{
			if (r<precision) 
				break;
			count=0;
		}
	}
	while (iter < nbSteps);

	return r;
	}


	//THE STUPID GSL CRASHES FROM TIME TO TIME!!!!!
	//BUT I USE IT CORRECTLY!

	ImplicitOnLine::start=start;
	ImplicitOnLine::diff= stop - start;
	ImplicitOnLine::imp=implicit;

	const gsl_root_fsolver_type *T=gsl_root_fsolver_bisection;
	gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);

	gsl_function F;

	F.function = ImplicitOnLine::proc;
	F.params =0;

	if (!s) return 0;
	double x_lo=0;
	double x_hi=1.0;
	double r;
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	gmVector3 recPos;
	unsigned int iter = 0;
	unsigned int count=0;
	int status;
	do
	{
		//this seems safer
		if (x_hi-x_lo<gmEPSILON)
			break;
		++iter;
		++count;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);

		//this seems a little too optimistic and often stops far from the surface
		//but even with it... the function crashes!
		if (status == GSL_SUCCESS)
			break;


		if (count==stepsPerPrecisionTest)
		{
			recPos=ImplicitOnLine::start+ ImplicitOnLine::diff*r;
			if (fabs(ImplicitOnLine::imp->proc(recPos))<precision) 
				break;
			count=0;
		}
	}
	while (iter < nbSteps);

	gsl_root_fsolver_free(s);
	return r;
}

void FindFirstIntersectionWithImplicits::drawParticle(int i)
{
	if (_shootTestRaysThroughParticles)
	{
		std::vector<Implicit * > implicits;
		implicits.push_back(_impInt->getImplicit());
		std::pair<gmVector3,Implicit*> result;
		findIntersection(	*_view->getCameraPosition(),_position->getPosition(i), implicits, &result, 
							_stepsPerUnit, _recursive, _precision, _stepsPerPrecisionTest);
		
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glDisable(GL_LIGHTING);
		glPointSize(10);
		glBegin(GL_POINTS);
		if (result.second)
		{
			glColor3f(0.5,1,0.5);
			glVertex3d(result.first[0],result.first[1],result.first[2]); 
		}
		else 
		{
			glColor3f(1,0.5,0.5);
			glVertex3d(_position->getPosition(i)[0],_position->getPosition(i)[1],_position->getPosition(i)[2]);
		}
		glEnd();
		glPopAttrib();
	}
}
