#include "ImplicitInterrogatorCached.h"
#include "ParticlePosition.h"
// why is Algebraic special?
#include "Surface/Implicit/Algebraic/Algebraic.h"
#include "ParticleSystem.h"		// for ->surfaces in attachAttributes()

REGISTER_PARTICLESTUFF(ImplicitInterrogatorCached,"Attribute:ImplicitInterrogatorCached");

/* Store name of implicit in qname[0]. Interface should present a menu of available
 * implicits indexed by their order, and accept one of these indices as setq(0).
 *
 * \todo If this routine is called with an implicit, its PSParam will reflect incorrectly that
 * the implicit is empty. -jch
 */

ImplicitInterrogatorCached::ImplicitInterrogatorCached(Particles* ps, Implicit* imp, const std::string &name)
	: ImplicitInterrogator(ps, imp, name)
{
	new PSParamBool(this,&_calculateProc,true,"proc","calc. proc","calculate the value of the implicit function at each particle");
	new PSParamDoublePerParticle(this,&_proc,"value","function value",
		"Value of the implicit surface function at this particle.");

	new PSParamBool(this,&_calculateGrad,true,"grad","calc. grad","calculate the gradient of the implicit function at each particle");
	new PSParamgmVector3PerParticle(this,&_grad,"grad","gradient",
		"Value of the gradient of the implicit function");

	new PSParamBool(this,&_calculateKg,true,"kg","calc. kg","calculate the gaussian curvature of the implicit function at each particle");
	new PSParamDoublePerParticle(this,&_kg,"kg","gaussian curvature",
		"Gaussian curvature measured at particle.");

	new PSParamBool(this,&_calculateKm,true,"km","calc. km","calculate the mean curvature of the implicit function at each particle");
	new PSParamDoublePerParticle(this,&_km,"km","mean curvature","Mean curvature measured at particle.");

	new PSParamBool(this, &_calculateHessian, true, "Hess","calc. Hess","Calculate the hessian of the implicit");	
	new PSParamDoublePerParticle(this,&_fxx,"Fxx","Fxx","Second partial derivative Fxx");
	new PSParamDoublePerParticle(this,&_fxy,"Fxy","Fxy","Second partial derivative Fxy");
	new PSParamDoublePerParticle(this,&_fxz,"Fxz","Fxz","Second partial derivative Fxz");
	new PSParamDoublePerParticle(this,&_fyy,"Fyy","Fyy","Second partial derivative Fyy");
	new PSParamDoublePerParticle(this,&_fyz,"Fyz","Fyz","Second partial derivative Fyz");
	new PSParamDoublePerParticle(this,&_fzz,"Fzz","Fzz","Second partial derivative Fzz");

}


/** Sets the particle system for which the
 * ImplicitInterrogator applies.
 * If new_ps is not NULL, then the lengths of internal arrays
 * of particle data are set to the size of the the particle
 * system.
**/
 void ImplicitInterrogatorCached::setParticleSystem(Particles *new_ps)
{
	ParticleAttribute::setParticleSystem(new_ps);
	if (ps)
	{
		_proc.resize(ps->size(), 0.0);
		_kg.resize(ps->size(), 0.0);
		_km.resize(ps->size(), 0.0);
		_grad.resize(ps->size(),gmVector3(0,0,0));
		_fxx.resize(ps->size(), 0.0);
		_fxy.resize(ps->size(), 0.0);
		_fxz.resize(ps->size(), 0.0);
		_fyy.resize(ps->size(), 0.0);
		_fyz.resize(ps->size(), 0.0);
		_fzz.resize(ps->size(), 0.0);
	}
}

/** Load caches
 */
void ImplicitInterrogatorCached::prepare()
{
	if (!_position) return;
	if (!_implicit) return;

	unsigned int i;
	unsigned int size=ps->size();

	if (_calculateProc)
	{
		for (i = 0; i < size; ++i)
			_proc[i] = _implicit->proc(_position->getPosition(i));
	}

	if (_calculateGrad)
	{
		for (i = 0; i < size; ++i)
			_grad[i] = _implicit->grad(_position->getPosition(i));
	}


	if (_calculateHessian)
	{
		for (i=0; i< size; ++i)
		{
			gmVector3 pos=_position->getPosition(i);
			gmMatrix3 h=_implicit->hess(pos);
	
			_fxx[i] = h[0][0];
			_fxy[i] = (h[0][1]+h[1][0])/2.0;
			_fxz[i] = (h[0][2]+h[2][0])/2.0;
			_fyy[i] = h[1][1];
			_fyz[i] = (h[1][2]+h[2][1])/2.0;
			_fzz[i] = h[2][2];
		}
	}

	if (_calculateKg)
	{
		//I don't think that there is any implicit not using Spivaks formula...
		if (_calculateGrad&&_calculateHessian)
		{
			for (i = 0; i < size; ++i)
			{
				gmVector3 gradient=grad(i);
				double fx=gradient[0];
				double fy=gradient[1];
				double fz=gradient[2];

				double den = gradient.lengthSquared();
				if (den == 0.0) 
					_kg[i]=HUGE;	// HUGE defined in math.h
				else
				{
					den = 1.0/den;

					double fxx = Fxx(i), fxy = Fxy(i), fxz = Fxz(i),
										 fyy = Fyy(i), fyz = Fyz(i),
													   fzz = Fzz(i);

					/* this is the version from sethian's level set book, but should probably
					* have |grad f|^4 in the denominator.
					*/
					/*
						double sethian_kg = ((fyy*fzz - fyz*fyz)*fx*fx + (fzz*fxx - fxz*fxz)*fy*fy + (fxx*fyy - fxy*fxy)*fz*fz +
							2.0*fx*fy*(fxz*fyz - fxy*fzz) - 2.0*fy*fz*(fxy*fxz - fyz*fxx) - 2.0*fx*fz*(fxy*fyz - fxz*fyy))
							* den;
					*/
					/* This is Spivak's version, taken from Suffern & Balsys 2003
						A = 2[fx fy(fxy fzz - fxz fyz) + fx fz(fxz fyy - fxy fyz) + fy fz(fyz fxx - fxy fxz)]
							- fx^2(fyy fzz - fyz^2) - fy^2(fxx fzz - fxz^2) - fz^2(fxxfyy - fxy^2)
						K = A/(fx^2 - fy^2 + fz^2)^2

						which should be -|hess(x) grad(x)|
										|grad(x)   0    |/|grad(x)|^4
					*/

					_kg[i]	=		(2.0*(fx*fy*(fxy*fzz - fxz*fyz) + fx*fz*(fxz*fyy - fxy*fyz) + fy*fz*(fyz*fxx - fxy*fxz)) -
									fx*fx*(fyy*fzz - fyz*fyz) - fy*fy*(fxx*fzz - fxz*fxz) - fz*fz*(fxx*fyy - fxy*fxy))
									* -den;
									//-(gradmag2 * gradmag2);
				}
			}
		}
		else
		{	
			for (i = 0; i < size; ++i)
				_kg[i] = _implicit->gaussianCurvature(_position->getPosition(i));
		}
	}

	if (_calculateKm)
	{
		if (_calculateGrad&&_calculateHessian)
		{
			for (i = 0; i < size; ++i)
			{
				gmVector3 gradient=grad(i);
				double fx = gradient[0];
				double fy = gradient[1];
				double fz = gradient[2];

				double fxx = Fxx(i), fxy = Fxy(i), fxz = Fxz(i),
									 fyy = Fyy(i), fyz = Fyz(i),
												   fzz = Fzz(i);

				_km[i] = 0.5*((fyy+fzz)*fx*fx + (fzz+fxx)*fy*fy + (fxx+fyy)*fz*fz -
						 2.0*fx*fy*fxy - 2.0*fy*fz*fyz - 2.0*fz*fx*fxz)/pow(fx*fx + fy*fy + fz*fz,1.5);
			}
		}
		else
		{
			for (i = 0; i < size; ++i)
				_km[i] = _implicit->meanCurvature(_position->getPosition(i));
		}
	}
}

/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see ParticleAttribute::particleRemoved
 *
 * \todo The default should manage per-particle array sizes automatically
 */
void ImplicitInterrogatorCached::particleRemoved(unsigned int i)
{

	// Copy end element to ith position
	_proc[i]  = _proc.back();
	_kg[i] = _kg.back();
	_km[i]  = _km.back();
	_grad[i]= _grad.back();
	_fxx[i]=_fxx.back();
	_fxy[i]=_fxy.back();
	_fxz[i]=_fxz.back();
	_fyy[i]=_fyy.back();
	_fyz[i]=_fyz.back();
	_fzz[i]=_fzz.back();

	// Pop off the end
	_proc.pop_back();
	_kg.pop_back();
	_km.pop_back();
	_grad.pop_back();
	_fxx.pop_back();
	_fxy.pop_back();
	_fxz.pop_back();
	_fyy.pop_back();
	_fyz.pop_back();
	_fzz.pop_back();
}

/** 
 * Callback for particle addition.
 * @param i Index of the particle that has been added.
 */
void ImplicitInterrogatorCached::particleAdded() 
{
	//maybe we want to calculate the correct new values?!
	//but for the moment we just leave it with one iteration of nonsense
	_proc.push_back(0.0);
	_kg.push_back(0.0);
	_km.push_back(0.0);
	_grad.push_back(gmVector3(0,0,0));
	_fxx.push_back(0.0);
	_fxy.push_back(0.0);
	_fxz.push_back(0.0);
	_fyy.push_back(0.0);
	_fyz.push_back(0.0);
	_fzz.push_back(0.0);
} 



double ImplicitInterrogatorCached::proc(unsigned int i) const
{
	return _proc[i];
}

gmVector3 ImplicitInterrogatorCached::grad(unsigned int i) const
{
	return _grad[i];
}


double ImplicitInterrogatorCached::Fx(unsigned int i) const
{
	return (grad(i))[0];
}
double ImplicitInterrogatorCached::Fy(unsigned int i) const
{
	return (grad(i))[1];
}
double ImplicitInterrogatorCached::Fz(unsigned int i) const
{
	return (grad(i))[2];
}

double ImplicitInterrogatorCached::Fxx(unsigned int i) const
{
	return _fxx[i];
}
double ImplicitInterrogatorCached::Fxy(unsigned int i) const
{
	return _fxy[i];
}
double ImplicitInterrogatorCached::Fxz(unsigned int i) const
{
	return _fxz[i];
}
double ImplicitInterrogatorCached::Fyy(unsigned int i) const
{
	return _fyy[i];
}
double ImplicitInterrogatorCached::Fyz(unsigned int i) const
{
	return _fyz[i];
}
double ImplicitInterrogatorCached::Fzz(unsigned int i) const
{
	return _fzz[i];
}

gmMatrix3 ImplicitInterrogatorCached::hess(unsigned int i) const
{
	//we exploit the symmetry of the hessian, 
	//if this causes numerical problems, we have to rethink this.
	//but in general I guess that having a symmetrical matrix, where symmetry
	//should be, will rather help.
	return gmMatrix3(_fxx[i], _fxy[i], _fxz[i],
					 _fxy[i], _fyy[i], _fyz[i],
					 _fxz[i], _fyz[i], _fzz[i]);
}


double ImplicitInterrogatorCached::kg(unsigned int i) const
{
	return _kg[i];
}

double ImplicitInterrogatorCached::km(unsigned int i) const
{
	return _km[i];
}
