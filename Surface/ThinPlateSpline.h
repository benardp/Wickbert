/**
* Declaration of an implicit 
 * @file ThinPlateSpline.h
 * @date 10 April. 2006
 * @author Eric Lorimer
 * @remarks
 */

#ifndef ThinPlateSpline_h
#define ThinPlateSpline_h

#include "Surface/Implicit/Implicit.h"
class ParticlePosition;
class ParticleScalar;


class ThinPlateSpline : public Implicit
{
public:
	ThinPlateSpline();	
	double proc(const gmVector3 & x);
	double proci(const gmVector3 & x, int i);
	gmVector3 grad(const gmVector3 & x);
	gmVector3 gradi(const gmVector3 & x, int i);
	
	~ThinPlateSpline() { }; 
	
    virtual void getq(double*);
    virtual void _setq(double*);
    virtual void getqname(char** qn);
	virtual unsigned int qlen() { return centers.size() * 4; }

	virtual void interpolate(Particles *ps, float phi);
	void updateRBF();
	double phi(const gmVector3 & x);

	MAKE_NAME();

protected:
	ParticlePosition *positions;
	ParticleScalar *values;

	std::vector<gmVector3> centers;
	std::vector<double> weights;
	double m_p[3];
};

#endif

