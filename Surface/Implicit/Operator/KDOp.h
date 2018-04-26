#ifndef		_KDOP_H
#define		_KDOP_H

#include "Surface/Implicit/Operator/UnaryOp.h"
#include "Surface/Implicit/Operator/KDTreeXin.h"
#include "Particles/Attributes/ParticlePosition.h"

class KDOpComp
{
public:
	double operator () (const gmVector3& pt, int i) {
		return pt[i];
	}
};

class KDOpTraits
{
public:
	typedef		gmVector3		point_type;
	typedef		double			coord_type;
	typedef		KDOpComp		get_component;
	static		const int		dim = 3;
};

/** KDOp is an operator that returns the distance from query point x to the closest particle position
 *  stored in the KD tree. If an implicit operand exists, returns the distance
 *  projected along the gradient of the operand implicit. KDOp is designed to accelerate the
 *  propagation of an implicit surface toward a set of points (stored as particle positions) by
 *  returning an Implicit used as the speed function and using the current procedural level set
 *  implicit surface approximation for the gradient to project the distance so that tangengtial
 *  distances result in zero speed.
 */
class KDOp : public UnaryOp
{
public:
	//	constructor
	KDOp();
	MAKE_NAME();

	//	prepare the surface
	void prepare();

	//	only overload the proc function
	virtual double proc(const gmVector3 & x);

protected:
	ParticlePosition*					position;
	//	kdtree
	xntools::ds::KDTree<KDOpTraits>		kd_tree;
};

/** KDOpParRefParam inherits all of the functionality of SurfAttrRefParam but adds
 * a call to KDOp::prepare() in KDOpAttrRefParam::attach(). This way the KDTree is
 * recomputed whenever the pointer to the particles representing its points changes.
 */
class KDOpAttrRefParam : public SurfAttrRefParam
{
public:
	KDOp *kdop;		///< Convenient cache of *parent, but of the right type

	KDOpAttrRefParam(KDOp *parent, ParticleAttribute **attr_ref, std::string s,
					std::string sn, std::string n, std::string d)
	: SurfAttrRefParam(parent,attr_ref, s, sn, n, d) { kdop = parent; }

	void attach(Surfaces *surfs, ParticleSystems *psystems) {
		SurfAttrRefParam::attach(surfs,psystems);
		kdop->prepare();
	}
};

#endif