/**
 * @file ABlend.h
 * @author John C. Hart
 * @date 15 Nov. 2005
 */

#ifndef ABLEND_H
#define ABLEND_H

#include "Blend.h"

/** ABlend is a pure virtual class used to collect shared functionality for the algebraic blends
 * described in Hoffmann & Hopcroft, Automatic Surface Generation in Computer Aided Design,
 * The Visual Computer 1, 1985, 92-100.
 * http://www.cs.purdue.edu/homes/cmh/distribution/papers/Geometry/geo1.pdf
 */
class ABlend : public Blend
{
  public:
    int m_d;    ///< Degree of blend

	/** if true, show just the blend surface (instead of the portion of the blend surface
	 * used to blend the two surfaces.
	 */
	bool blendsurf;

    /// Constuctor defaults to a quadratic with radii one.
	ABlend(Implicit *f = NULL, Implicit *g = NULL, double r1 = 1.0, double r2 = 1.0, int d = 2) :
		Blend(f,g,r1,r2)
	{
		new SurfParamInt(this,&m_d,d,"deg","degree","Degree of blend.");
		new SurfParamBool(this,&blendsurf,false,"blendsurf","Blend Surface","Show the blend surface only.");
	}

	/// Blending region is limited to f < r1, g < r2
	virtual bool indomain(double f, double g) {
		return blendsurf ? true : (fabs(f) <= m_r1 && fabs(g) <= m_r2);
	}

	/** 
	 * Computes an elliptic blend.
	 * h(f,g) = -( (f - r1)^d/r1^d + (g - r2)^d/r2^d - 1).
	 * Negated such that surface is negative inside.
	 * Returns smaller of f or g if either is outside its
	 * blending radius.
	 *
	 * When f = 0 and as g approaches r2, returns zero.
	 * When g = 0 and as f approaches r1, returns zero.
	 *
	 * To get intersection, compute the complement of the union of the complement
	 */
	virtual double h(double f, double g) {
		return indomain(f,g) ? 1.0 - pow(fabs(1.0-f/m_r1),m_d) - pow(fabs(1.0-g/m_r2),m_d) : (f<g) ? f : g;
	}

	///dh(f,g)/df = (- -) d (1-f/r1)^(d-1)
	virtual double hf(double f, double g) {
		return indomain(f,g) ? m_d*pow(1.0 - f/m_r1,m_d-1.0)/m_r1 : (f < g) ? 1.0 : 0.0;
	}

	/// dh(f,g)/dg = (- -) d (1-g/r2)^(d-1)
	virtual double hg(double f, double g) {
		return indomain(f,g) ? m_d*pow(1.0 - g/m_r2,m_d-1.0)/m_r2 : (f < g) ? 0.0 : 1.0;
	}

	/// d2h(f,g)/df2 = - d(d-1) (1-f/r1)^(d-2), (m_d should be at least two?)
	virtual double hff(double f, double g) {
		return indomain(f,g) ? -m_d*(m_d-1.0)*pow(1.0 - f/m_r1,m_d-2.0)/(m_r1*m_r1) : 0.0;
	}

	/// d2h(f,g)/dg2 = - d(d-1) (1-g/r2)^(d-2)
	virtual double hgg(double f, double g) {
		return indomain(f,g) ? -m_d*(m_d-1.0)*pow(1.0 - g/m_r2,m_d-2.0)/(m_r2*m_r2) : 0.0;
	}

	/// No mixed second derivatives
	double hfg(double f, double g) { return 0.0; }


	/* ...INTERVAL VERSIONS... */


	/** Not sure if we really need this to be interval, but here it is.
	 ** should probably return a boolean interval, which is probably tristate
	 */
	virtual bool indomain(Intervald f, Intervald g) {
		return f.overlaps(Intervald(0.0,m_r1)) && g.overlaps(Intervald(0.0,m_r2));
	}

	/** Interval version of blending function */
	virtual Intervald h(Intervald f, Intervald g)
	{
		if (indomain(f,g))
			return Intervald(1.0) - 
				(Intervald(1.0) - f/Intervald(m_r1)).pow(Intervald(m_d)) - 
				(Intervald(1.0) - g/Intervald(m_r2)).pow(Intervald(m_d));
		else
			return (f < g) ? f : g;
	}

	virtual Intervald hf(Intervald f, Intervald g)
	{
		if (indomain(f,g))
			return Intervald(m_d) * 
				(Intervald(1.0) - f/Intervald(m_r1)).pow(Intervald(m_d-1.0)) / Intervald(m_r1);
		else
			return (f < g) ? Intervald(1.0) : Intervald(0.0);
	}

	virtual Intervald hg(Intervald f, Intervald g)
	{
	if (indomain(f,g))
		return Intervald(m_d) *
			(Intervald(1.0) - g/Intervald(m_r2)).pow(Intervald(m_d-1.0)) / 
		Intervald(m_r2);
	else
		return (f < g) ? Intervald(0.0) : Intervald(1.0);
	}

	virtual Intervald hff(Intervald f, Intervald g)
	{
		if (indomain(f,g))
			return Intervald(m_d) * Intervald(m_d-1.0) * 
				(Intervald(1.0) - f/Intervald(m_r1)).pow(Intervald(m_d-2.0)) / Intervald(m_r1*m_r1);
		else
			return Intervald(0.0);
	}

	virtual Intervald hgg(Intervald f, Intervald g)
	{
	if (indomain(f,g))
		return Intervald(m_d) * Intervald(m_d-1.0) *
			(Intervald(1.0) - g/Intervald(m_r2)).pow(Intervald(m_d-2.0)) / 
			Intervald(m_r2*m_r2);
	else
		return Intervald(0.0);
	}

	Intervald hfg(Intervald f, Intervald g) { return Intervald(0.0); }
};

/** AUnion just implements ABlend which defaults to a blended union
 */
class AUnion : public ABlend
{
public:
    /// Constuctor defaults to a quadratic with radii one.
	AUnion(Implicit *f = NULL, Implicit *g = NULL, double r1 = 1.0, double r2 = 1.0, int d = 2) :
	  ABlend(f,g,r1,r2,d) { }

    MAKE_NAME();
};

/** AIntersection implements intersection as the complement of the union of the complements.
 * Note that the derivatives are wrt f and g which when negated bring an extra (-) sign out
 * with each derivative.
 */
class AIntersection : public ABlend
{
public:
    /// Constuctor defaults to a quadratic with radii one.
	AIntersection(Implicit *f = NULL, Implicit *g = NULL, double r1 = 1.0, double r2 = 1.0, int d = 2) :
	  ABlend(f,g,r1,r2,d) { }

	/// intersection is complement of the union of the complements
	virtual double h(double f, double g) { return -ABlend::h(-f,-g); }
	virtual double hf(double f, double g) { return ABlend::hf(-f,-g); }
	virtual double hg(double f, double g) { return ABlend::hg(-f,-g); }
	virtual double hff(double f, double g) { return -ABlend::hff(-f,-g); }
	virtual double hgg(double f, double g) { return -ABlend::hgg(-f,-g); }

    MAKE_NAME();
};

/** ASubtraction implements subtraction (intersection of A with the complement of B) as
 * -(-A U -(-B)) = -(-A U B).
 */
class ASubtraction : public ABlend
{
public:
    /// Constuctor defaults to a quadratic with radii one.
	ASubtraction(Implicit *f = NULL, Implicit *g = NULL, double r1 = 1.0, double r2 = 1.0, int d = 2) :
	  ABlend(f,g,r1,r2,d) { }

	/// intersection is complement of the union of the complements
	virtual double h(double f, double g) { return -ABlend::h(-f,g); }
	virtual double hf(double f, double g) { return ABlend::hf(-f,g); }
	virtual double hg(double f, double g) { return -ABlend::hg(-f,g); }
	virtual double hff(double f, double g) { return -ABlend::hff(-f,g); }
	virtual double hgg(double f, double g) { return -ABlend::hgg(-f,g); }

    MAKE_NAME();
};

#endif

