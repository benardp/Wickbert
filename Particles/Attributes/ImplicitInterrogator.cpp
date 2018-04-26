#include "ImplicitInterrogator.h"
#include "ParticlePosition.h"
// why is Algebraic special?
#include "Surface/Implicit/Algebraic/Torus.h"
#include "ParticleSystem.h"		// for ->surfaces in attachAttributes()

REGISTER_PARTICLESTUFF(ImplicitInterrogator,"Attribute:ImplicitInterrogator");

/** empty_imp is a default implicit that the implicit interrogator points to when
 * there are no other implicits at which to point.
 * Pointed to by the global Implicit* empty.
 */
/** empty is a pointer to a default implicit f(x)=0.
 * Points to the Algebraic empty_imp.
 */
static Implicit *empty_p = dynamic_cast<Implicit*>(new Torus());
//Algebraic(0)


/* Store name of implicit in qname[0]. Interface should present a menu of available
 * implicits indexed by their order, and accept one of these indices as setq(0).
 *
 * \todo If this routine is called with an implicit, its PSParam will reflect incorrectly that
 * the implicit is empty. -jch
 */

ImplicitInterrogator::ImplicitInterrogator(Particles* ps, Implicit* imp, const std::string &name)
	: ParticleAttribute(ps, name)
{
	new PSParamImplicitName(this,&impname,"<empty>",
		"implicit","implicit surface","Name of the implicit surface");

	new Attached<ParticlePosition>(this, &_position);

	index=0;
	if (imp)
		_implicit=imp;
	else
		empty();
}

void ImplicitInterrogator::attachAttributes()
{
	ParticleAttribute::attachAttributes();
	setImplicit(impname);
}

/** Sets the particle system for which the
 * ImplicitInterrogator applies.
 * If new_ps is not NULL, then the lengths of internal arrays
 * of particle data are set to the size of the the particle
 * system.
 *
 * \todo The default setParticleSystem should resize per-particle arrays automatically
 */
void ImplicitInterrogator::setParticleSystem(Particles *new_ps)
{
	ParticleAttribute::setParticleSystem(new_ps);
}


void ImplicitInterrogator::setImplicit(Implicit* imp)
{ 
	//assert(implicit);
	_implicit = imp;
	if (imp)
	impname = imp->getObjectName();
	//if (implicit)
	//	dqdt.resize(0, implicit->qlen());
}

void ImplicitInterrogator::setImplicit(std::string s)
{
	_implicit = empty_p;
	ParticleSystem *psys = ps->particleSystem;
	if (psys) {
		Surfaces *surfs = psys->surfaces;
		if (surfs) {
			_implicit = dynamic_cast<Implicit *>(surfs->withName(s));
			if (!_implicit) _implicit = empty_p;
		}
	}
}

void ImplicitInterrogator::empty()
{
	setImplicit(empty_p);
}


double ImplicitInterrogator::proc(unsigned int i) const
{
	return _implicit->proc(_position->getPosition(i));
}

gmVector3 ImplicitInterrogator::grad(unsigned int i) const
{
	return _implicit->grad(_position->getPosition(i));
}

gmVector3 ImplicitInterrogator::normal(unsigned int i)const
{
	//copied from implicit
	gmVector3 n = grad(i);
  
	if (gmIsZero(n.length()))
		n = gmVector3(1.0,0.0,0.0);
	else
		n.normalize();
  
	return n;
}

double ImplicitInterrogator::Fx(unsigned int i) const
{
	return _implicit->Fx(_position->getPosition(i));
}
double ImplicitInterrogator::Fy(unsigned int i) const
{
	return _implicit->Fy(_position->getPosition(i));
}
double ImplicitInterrogator::Fz(unsigned int i) const
{
	return _implicit->Fz(_position->getPosition(i));
}
double ImplicitInterrogator::Fxx(unsigned int i) const
{
	return _implicit->Fxx(_position->getPosition(i));
}
double ImplicitInterrogator::Fxy(unsigned int i) const
{
	return _implicit->Fxy(_position->getPosition(i));
}
double ImplicitInterrogator::Fxz(unsigned int i) const
{
	return _implicit->Fxz(_position->getPosition(i));
}
double ImplicitInterrogator::Fyy(unsigned int i) const
{
	return _implicit->Fyy(_position->getPosition(i));
}
double ImplicitInterrogator::Fyz(unsigned int i) const
{
	return _implicit->Fyz(_position->getPosition(i));
}
double ImplicitInterrogator::Fzz(unsigned int i) const
{
	return _implicit->Fzz(_position->getPosition(i));
}
gmMatrix3 ImplicitInterrogator::hess(unsigned int i) const
{
	return _implicit->hess(_position->getPosition(i));
}

double ImplicitInterrogator::kg(unsigned int i) const
{
	return _implicit->gaussianCurvature(_position->getPosition(i));
}

double ImplicitInterrogator::km(unsigned int i) const
{
	return _implicit->meanCurvature(_position->getPosition(i));
}
