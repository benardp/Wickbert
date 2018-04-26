/**
 * Implementation of the particle visibility attribute.
 * @file ParticleVisibility.cpp
 * @author Matei N. Stroila
 * @date 6/10/2005
 * @remarks 
 */

#include "ParticleVisibility.h"

#include "ParticleBoundingBox.h"
#include "ParticlePosition.h"
#include "ParticleNormal.h"
#include "ViewDependence.h"
#include "ImplicitInterrogator.h"

#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"

#include <algorithm>


REGISTER_PARTICLESTUFF(ParticleVisibility,"Attribute:ParticleVisibility");


ParticleVisibility::VisibilityMethodCallback::VisibilityMethodCallback(ParticleVisibility * pv)
	:PSParamComboBox::Callback()
	,_pv(pv)
{
	_choices.push_back(std::string("GSL"));
	_choices.push_back(std::string("STEPPING"));
}

void ParticleVisibility::VisibilityMethodCallback::itemselected(unsigned int selection)
{
	if (selection==0)//GSL
		_pv->visibilityMethod=ParticleVisibility_GSL;
	else
		_pv->visibilityMethod=ParticleVisibility_STEPPING;
}


/**
* Adds ParticleVisibility to a system of particles.
 * \param ps   The owning particle system.
 * \param name The name of this object.
 */
ParticleVisibility::ParticleVisibility(Particles *ps, const std::string& name)
: ParticleAttribute(ps, name) 
, visibilityMethod(ParticleVisibility_GSL)
{
	if (ps)
		mVisibility.resize(ps->size(),1);
	new Attached<ParticlePosition>(this,&pos);
	new Attached<ParticleNormal>(this,&orient);
	new Attached<ViewDependence>(this,&viewDep);
	new Attached<ParticleBoundingBox>(this,&bbox);

	new PSParamBool(this,&deactivated,false,"deact","deactivate","No visibility calculation is performed. All particles are visible.");
	isDeactivated=false;
	new PSParamBool(this,&lazyEvaluation,false,"lazy","lazy Ev.:","Lazy evaluation, only if the camera moves.");
	new PSParamBool(this,&useNormalTest,false,"useNormalTest","Test normals","Set true if want to use backface culling, should be false for silhouettes");
	new PSParamBool(this,&useBbox,false,"useBbox","Use bounding box","Set true if want to use bounding box intersections");
	new PSParamDouble(this,&bboxOffset,1.5,"bboxOffset","Bbox offset","Bbox offset");
	new PSParamComboBox(this,new VisibilityMethodCallback(this), "visMethod", "visibilityMethod", "Choose the method used for visibility test");

	new PSParamDouble(this,&newtonEps,0.1,"newtonEps","Newton Epsilon:GSL", "Epsilon for the Newton's method");
	new PSParamDouble(this,&offset,0.0,"offset","Offset:GSL","Offset");
	new PSParamDouble(this,&initX,0.9,"initX","initX:GSL" ,"Where to init Newton iterations (0.0 = particle, 1.0 = eye)");
	new PSParamDouble(this,&minX,0.01,"minX","minX:GSL" ,"Threshold for the Newton solver");
	new PSParamInt(this,&maxIter,10,"maxIter","Max. no. of Iterations:GSL", "GSLParam: Maximum number of Iterations");
	new PSParamDouble(this,&nbStepsPerUnit,32.0,"nbStepsPerUnit","nbStepsPerUnit:STEP","Stepping Param: Nr steps per unit");

	new PSParamIntPerParticle(this,&mVisibility,"vis","visibility", "1 indicates this particle is not occluded from view");

} // end ParticleVisibility::ParticleVisibility()


void ParticleVisibility::prepare()
{
	//I don't think we need this...
	//if (!ps) return;
	//if (!(ps->particleSystem)) return;
	//if (!(ps->particleSystem->particleSystems))return;
	assert(ps);
	assert(ps->particleSystem);
	assert(ps->particleSystem->particleSystems);

	//static unsigned int t=0;
	//++t;
	//if (t==100)
	//{
	//	if (useBbox==true)
	//		std::cout<<"BOUNDINGBOX"<<'\n';
	//	else if (useBbox==false)
	//		std::cout<<"no bounding!"<<'\n';
	//	if (lazyEvaluation==true)
	//		std::cout<<"LAZY"<<'\n';
	//	else if (lazyEvaluation==true)
	//		std::cout<<"no lazy"<<'\n';
	//		
	//	if (deactivated==true)
	//		std::cout<<"DEACTIVATED!"<<'\n';
	//	else if (deactivated==false)
	//		std::cout<<"activated"<<'\n';
	//	if (visibilityMethod==ParticleVisibility_GSL)
	//		std::cout<<"GSL"<<'\n';
	//	else if (visibilityMethod==ParticleVisibility_STEPPING)
	//		std::cout<<"STEPPING"<<'\n';
	//	else
	//		std::cout<<"DEEP TROUBLE!!!"<<'\n';
	//	t=0;
	//}



	if (deactivated)
	{
		//if the set is not already reset, then let's do it.
		if (!isDeactivated)
		{
			for (unsigned int i=0;i<mVisibility.size();++i)
				mVisibility[i]=1;
			
			isDeactivated=true;
		}
		return;
	}

	isDeactivated=false;

	//if lazy we only reevaluate if camera moves a lot.
	if (lazyEvaluation)
	{
		if ((*(viewDep->getCameraPosition())-oldCameraPosition).length()<gmEPSILON)
			return;
		oldCameraPosition=*(viewDep->getCameraPosition());
	}
	
	_implicits.clear();
	

	
	//for each implicit in the opengl scene, check the visibility
	std::vector<Particles*>::iterator piter;
	ParticleSystems::iterator pSysIter;

	for(pSysIter = (*(ps->particleSystem->particleSystems)).begin(); pSysIter != (*(ps->particleSystem->particleSystems)).end();  ++pSysIter)
	{
		for(piter = ((*pSysIter)->particles).begin(); piter != ((*pSysIter)->particles).end(); ++piter)
		{	
			std::vector<ImplicitInterrogator *> interrogators;

			(*piter)->getAllAttributesByType<ImplicitInterrogator>(&interrogators);
			for (unsigned int i=0;i<interrogators.size();++i)
			{
				Implicit* imp=interrogators[i]->getImplicit();
				if (imp!=0)
					if (find(_implicits.begin(),_implicits.end(),imp)==_implicits.end())
						_implicits.push_back(imp);
			}
		}
	}
	
	//std::cout<<"found the following"<<_implicits.size()<<"implicits"<<std::endl;
	//for (unsigned int i=0;i<_implicits.size();++i)
	//	std::cout<<_implicits[i]->name()<<std::endl;
	
	//_implicits.clear();

	//for (Surfaces::iterator iter=ps->particleSystem->surfaces->begin();
	//	iter!=ps->particleSystem->surfaces->end(); ++iter)
	//{
	//	Implicit* imp;
	//	if (imp=dynamic_cast<Implicit*>(*iter))
	//		_implicits.push_back(imp);
	//}
	
	//std::cout<<"IN Program found the following"<<_implicits.size()<<"implicits"<<std::endl;
	//for (unsigned int i=0;i<_implicits.size();++i)
	//	std::cout<<_implicits[i]->name()<<std::endl;

	//std::cout<<"_________________________"<<std::endl;
	

	if (visibilityMethod==ParticleVisibility_GSL)
	{
		//we don't want the program to crash!
		gsl_set_error_handler_off();

		for (unsigned int i=0;i<ps->size();++i)
		{
			if(useNormalTest)
			{
				double _dot = dot(orient->getNormal(i),*(viewDep->getCameraPosition()) - pos->getPosition(i));
				if(_dot < -gmEPSILON)	
				{
					mVisibility[i]=0;
				}
				else
					mVisibility[i] = computeVisibilityGSL(i);
			}
			else
				mVisibility[i] = computeVisibilityGSL(i);
		}
	}
	else 
	{
		for (unsigned int i=0;i<ps->size();++i)
		{
			if(useNormalTest)
			{
				double _dot = dot(orient->getNormal(i),*(viewDep->getCameraPosition()) - pos->getPosition(i));
				if(_dot < -gmEPSILON)	
				{
					mVisibility[i]=0;
				}
				else
					mVisibility[i] = computeVisibilityStepping(i);
			}
			else
				mVisibility[i] = computeVisibilityStepping(i);
		}
	}
}



void ParticleVisibility::reset()
{
	for(unsigned int i=0; i<mVisibility.size(); i++)
		mVisibility[i] = 1;
} // end ParticleVisibility::reset()

void ParticleVisibility::clear()
{ 
	mVisibility.clear();
} // end ParticleVisibility::clear()

/**
* Add a visibility corresponding to the new particle.
 * \param i Index of the new particle.
 */
void ParticleVisibility::particleAdded() 
{
	mVisibility.push_back(1);
} // end ParticleVisibility::particleAdded()

/**
* Callback for particle removal.
 * \param i Index of particle to be removed.
 * \see Particles::particleRemoved
 */
void ParticleVisibility::particleRemoved(unsigned int i) 
{
	mVisibility[i] = mVisibility.back();
	mVisibility.pop_back();
} // end ParticleVisibility::particleRemoved()

int ParticleVisibility::getVisibility( gmVector3* cameraPosition,  unsigned int i)
{
	//debug...
	//because this should not be possible to call i with camera position...
	//Please remove all the calls of this function.
	//the lazy evaluation should be part of THIS class, not of the others. 
	//I put it in as an option.
	return mVisibility[i];	
	
} // end ParticleVisibility::getVisibility(const gmVector3& cameraPosition, const unsigned int i)


int ParticleVisibility::computeVisibilityStepping(int i) const
{	
	gmVector3 intersectBbox;
	if(useBbox)
		bbox->computeIntersection(i,*(viewDep->getCameraPosition()),
									bboxOffset,
									intersectBbox );
	else
		intersectBbox = *(viewDep->getCameraPosition());


	for (unsigned int currImp=0;currImp<_implicits.size();++currImp)
	{
			gmVector3 viewx(pos->getPosition(i) - intersectBbox );
			unsigned int nbSteps = (unsigned int) (viewx.length()* nbStepsPerUnit);
			for (unsigned int i=0;i<nbSteps;++i)
			{
				gmVector3 evalV=viewx * ((double)i)/((double)nbSteps) + intersectBbox;
				double val=_implicits[currImp]->proc(evalV);
				if (val<0)
					return 0;
			}
	}
	return 1;
}
int ParticleVisibility::computeVisibilityGSL(int i) const
{
	//offset for silhouettes
	gmVector3 offsetVector = *(viewDep->getCameraPosition()) - pos->getPosition(i);
	offsetVector =  pos->getPosition(i) + offset * offsetVector.normalize();
	
	gmVector3 intersectBbox;
	if(useBbox)
		bbox->computeIntersection(i,*(viewDep->getCameraPosition()),
									bboxOffset,
									intersectBbox );
	else
		intersectBbox = *(viewDep->getCameraPosition());


	
	for (unsigned int currImp=0;currImp<_implicits.size();++currImp)
	{
		struct myF_params params = {_implicits[currImp], &intersectBbox, &offsetVector, bbox};
		int status;
		int iter = 0;
		double x0, x = initX;
		
		const gsl_root_fdfsolver_type *T;
		gsl_root_fdfsolver *s;
		gsl_function_fdf FDF;
		FDF.f = &myF;
		FDF.df = &myF_deriv;
		FDF.fdf = &myF_fdf;
		FDF.params = &params;
		T = gsl_root_fdfsolver_steffenson;
		s = gsl_root_fdfsolver_alloc (T);	
		gsl_root_fdfsolver_set (s, &FDF, x);
		do
		{
			iter++;
			status = gsl_root_fdfsolver_iterate (s);
			x0 = x;
			x = gsl_root_fdfsolver_root (s);
			double _f = myF(x,(void*) &params);
			status =  gsl_root_test_residual (_f, newtonEps);
		}
		while (status == GSL_CONTINUE && iter < maxIter );
		
		gsl_root_fdfsolver_free(s);	
		if( x < initX && x > minX && status == GSL_SUCCESS )
			return 0;
	}
	return 1;
}// end ParticleVisibility::computeVisibility(int i)

//functions for GSL
//the viewVector = cameraPosition - particlePosition


//WHY DO WE NEED THE INBOUNDS TESTS????
// IT ALSO WORKS WITHOUT AND I SEE NO REASON TO ADD THEM...
//I comment them out.
//it would be nice to have it without them and I see no reason why not.
//then the ParticleVisibility could have one public const function:
//bool testSegmentForIntersection(s1,s2);
//it is basically a copy paste of the above.
//and then we call in the above exactly this function.
double myF (double x, void *params)
{
	struct myF_params *p = (struct myF_params *) params;
	gmVector3 viewx((*(p->myIntersectBbox) - *(p->myParticlePos)) * x +  *(p->myParticlePos));
//	if(p->bbox->inBounds(viewx))
		return p->myImp->proc(viewx);	
//	else
//		return gmGOOGOL;
}
double myF_deriv (double x, void *params)
{
	struct myF_params *p = (struct myF_params *) params;
	gmVector3 viewDirection(*(p->myIntersectBbox) - *(p->myParticlePos));
	gmVector3 viewx(viewDirection * x +  *(p->myParticlePos));	
//	if(p->bbox->inBounds(viewx))
		return dot(p->myImp->grad(viewx),viewDirection);
//	else
//		return gmGOOGOL;
}
void myF_fdf (double x, void *params, 
				 double *y, double *dy)
{
	struct myF_params *p  = (struct myF_params *) params;
	gmVector3 viewDirection(*(p->myIntersectBbox) - *(p->myParticlePos));
	gmVector3 viewx(viewDirection * x +  *(p->myParticlePos));	
//	if(p->bbox->inBounds(viewx)){
		*y = p->myImp->proc(viewx);
		*dy =  dot(p->myImp->grad(viewx),viewDirection);
//	}
//	else
	//{
	//	*y = gmGOOGOL;
	//	*dy = gmGOOGOL;
	//}
}
