
#include "KDOp.h"
#include <fstream>
#include <sstream>
using namespace std;

REGISTER_IMPLICIT(KDOp, "UnaryOp:KDOp");

//	construct the function with a vector of gmVector3
KDOp::KDOp()
{
	new KDOpAttrRefParam(this,(ParticleAttribute **)&position,"p","Points",
			"Particles attribute containing the data point positions", "Points");
}

void KDOp::prepare()
{
	//	construct the kd tree
	if (!position)
		return;
	
	kd_tree.construct(position->x);
}

/** Returns the distance from query point x to the closest particle position
 *  stored in the KD tree. If an implicit operand exists, returns the distance
 *  projected along the gradient of the operand implicit.
 */
double KDOp::proc(const gmVector3 & x)
{
	if (!position) return 0.0;

	double v[3];
	x.copyTo(v);
	int id = kd_tree.search_nearest(v, position->x);
	gmVector3 dist = position->x[id] - x;
	if (m_f) {
		gmVector3 grad = m_f->grad(x);
		grad.normalize();
		return dot(dist,grad);
	} else {
		return dist.length();
	}
}
