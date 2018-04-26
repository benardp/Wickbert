/**
 * @file RFunction.h
 * Common base class for unions and intersections; both are defined using
 * R-Functions, with only a difference in sign in one of the terms.
 *
 * @author Steve Zelinka
 */

#ifndef _RFUNCTION_H
#define _RFUNCTION_H

#include "Blend.h"

/** 
 * A subclass of blending that implements R-functions for CSG.
 * R-functions are blending functions that apply CSG operations
 * to implicit surfaces. The benefit of R-functions are that they
 * smooth the field surrounding the zero-surfaces to eliminate
 * the creasing that can result from, say, min and max functions.
 */
class RFunction : public BinaryOp
{
public:
	RFunction(double cont = 1.0, double sign = -1.0);
  protected:
    /**
     * The continuity to be provided.
     * This is stored as a double to avoid integer
     * round-offs later.
	 * 
	 * Should be a SurfParamInt -jh
     */
    double m_cont;

    /**
     * The sign of the second-half of the RFunction.
     * This differentiates
     * a union (= -1) from an intersection (= 1).
     * @note The sign is opposite of that used in the
     * original R-function paper where f<0 outside.
     */
    double m_sign;

    /**
     * Evaluates the R-function on the two function values.
     *
     * @param f  The value of f.
     * @param g  The value of g.
     * @return (f + g + m_sign*(sqrt(f^2+g^2)))*(f^2+g^2)^(m_cont/2).
     */
    virtual double h(double f, double g);
    virtual Intervald h(Intervald f, Intervald g);

    /**
     * Evaluates the partial derivative with respect to f
     * on the two function values.
     *
     * h = (f + g + m_sign*(sqrt(f^2+g^2)))*(f^2+g^2)^(m_cont/2).
     * hf = (f^2+g^2)^(m_cont/2)*d(f + g + m_sign*(sqrt(f^2+g^2)))/df +
     *      (f + g + m_sign*(sqrt(f^2+g^2)))*d(f^2+g^2)^(m_cont/2)/df
     *    = (f^2+g^2)^(m_cont/2)*(1 + m_sign*f/sqrt(f^2+g^2)) +
     *      (f + g + m_sign*(sqrt(f^2+g^2)))*m_cont*f*(f^2+g^2)^(m_cont/2 - 1).
     *
     * @param f The value of f.
     * @param g The value of g.
     * @return hf.
     *
     * \note Appears to work fine when tested with particles.
     */
    virtual double hf(double f, double g);
    virtual Intervald hf(Intervald f, Intervald g);

    /**
     * Evaluates the partial derivative with respect to g
     * on the two function values.
     *
     * h = (f + g + m_sign*(sqrt(f^2+g^2)))*(f^2+g^2)^(m_cont/2).
     * hg = (f^2+g^2)^(m_cont/2)*d(f + g + m_sign*(sqrt(f^2+g^2)))/dg +
     *      (f + g + m_sign*(sqrt(f^2+g^2)))*d(f^2+g^2)^(m_cont/2)/dg
     *    = (f^2+g^2)^(m_cont/2)*(1 + m_sign*g/sqrt(f^2+g^2)) +
     *      (f + g + m_sign*(sqrt(f^2+g^2)))*m_cont*g*(f^2+g^2)^(m_cont/2 - 1).
     *
     * @param f The value of f.
     * @param g The value of g.
     * @return hg.
     *
     * \note Appears to work fine when tested with particles.
     */
    virtual double hg(double f, double g);
    virtual Intervald hg(Intervald f, Intervald g);

    /**
     * Evaluates the 2nd partial derivative with respect to f
     * and f on the two function values.
     *
     * fg = sqrt(f^2 + g^2)
     * dfg/df = f/fg
     * hf  = fg^m_cont*(1 + m_sign*f/fg) +
     *       (f + g + m_sign*fg)*m_cont*f*fg^(m_cont - 2)
     *     = fg^m_cont + m_sign*f*fg^(m_cont-1) +
     *       (f + g + m_sign*fg)*m_cont*f*fg^(m_cont - 2)
     * hff = m_cont*fg^(m_cont-1)*dfg/df +
     *       m_sign*f*(m_cont-1)*fg^(m_cont-2)*dfg/df +
     *       m_sign*fg^(m_cont-1) +
     *       (1 + m_sign*dfg/df)*m_cont*f*fg^(m_cont - 2) +
     *       (f + g + m_sign*fg)*m_cont*(f*(m_cont - 2)*fg^(m_cont - 3)*dfg/df +
     *       fg^(m_cont - 2))
     *     = m_cont*f*fg^(m_cont-2) +
     *       m_sign*(m_cont-1)*f^2*fg^(m_cont-3) +
     *       m_sign*fg^(m_cont-1) +
     *       m_cont*f*fg^(m_cont - 2) + m_cont*m_sign*f^2*fg^(m_cont - 3) +
     *       m_cont*(f + g + m_sign*fg)*((m_cont - 2)*f^2*fg^(m_cont - 4) + 
     *       fg^(m_cont - 2))
     * 
     * @param f The value of f.
     * @param g The value of g.
     * @return hff.
     *
     * \note Needs to be tested.
     */
    virtual double hff(double f, double g);
    virtual Intervald hff(Intervald f, Intervald g);

    /**
     * Evaluate the 2nd partial derivative with respect to f
     * and g on the two function values.
     *
     * fg = sqrt(f^2 + g^2)
     * dfg/dg = g/fg
     * hf  = fg^m_cont*(1 + m_sign*f/fg) +
     *       (f + g + m_sign*fg)*m_cont*f*fg^(m_cont - 2)
     *     = fg^m_cont + m_sign*f*fg^(m_cont-1) +
     *       (f + g + m_sign*fg)*m_cont*f*fg^(m_cont - 2)
     * hfg = m_cont*fg^(m_cont-1)*dfg/dg +
     *       m_sign*f*(m_cont-1)*fg^(m_cont-2)*dfg/dg +
     *       m_sign*fg^(m_cont-1) +
     *       (1 + m_sign*dfg/dg)*m_cont*f*fg^(m_cont - 2) +
     *       (f + g + m_sign*fg)*m_cont*(f*(m_cont - 2)*fg^(m_cont - 3)*dfg/dg +
     *       fg^(m_cont - 2))
     *     = m_cont*g*fg^(m_cont-2) +
     *       m_sign*(m_cont-1)*f*g*fg^(m_cont-3) +
     *       m_sign*fg^(m_cont-1) +
     *       m_cont*f*fg^(m_cont - 2) + m_cont*m_sign*f*g*fg^(m_cont - 3) +
     *       m_cont*(f + g + m_sign*fg)*((m_cont - 2)*f*g*fg^(m_cont - 4) + 
     *       fg^(m_cont - 2))
     * 
     * @param f The value of f.
     * @param g The value of g.
     * @return hfg.
     *
     * \note Needs to be tested.
     */
    virtual double hfg(double f, double g);
    virtual Intervald hfg(Intervald f, Intervald g);

    /**
     * Return the 2nd partial derivative with respect to g
     * and g on the two function values.
     *
     * fg = sqrt(f^2 + g^2)
     * dfg/dg = g/fg
     * hg  = fg^m_cont*(1 + m_sign*g/fg) +
     *       (f + g + m_sign*fg)*m_cont*g*fg^(m_cont - 2)
     *     = fg^m_cont + m_sign*g*fg^(m_cont-1) +
     *       (f + g + m_sign*fg)*m_cont*g*fg^(m_cont - 2)
     * hgg = m_cont*fg^(m_cont-1)*dfg/dg +
     *       m_sign*g*(m_cont-1)*fg^(m_cont-2)*dfg/dg +
     *       m_sign*fg^(m_cont-1) +
     *       (1 + m_sign*dfg/dg)*m_cont*g*fg^(m_cont - 2) +
     *       (f + g + m_sign*fg)*m_cont*(g*(m_cont - 2)*fg^(m_cont - 3)*dfg/dg +
     *       fg^(m_cont - 2))
     *     = m_cont*g*fg^(m_cont-2) +
     *       m_sign*(m_cont-1)*g^2*fg^(m_cont-3) +
     *       m_sign*fg^(m_cont-1) +
     *       m_cont*g*fg^(m_cont - 2) + m_cont*m_sign*g^2*fg^(m_cont - 3) +
     *       m_cont*(f + g + m_sign*fg)*((m_cont - 2)*g^2*fg^(m_cont - 4) + 
     *       fg^(m_cont - 2))
     * 
     * @param f The value of f.
     * @param g The value of g.
     * @return hgg.
     *
     * \note Needs to be tested.
     */
    virtual double hgg(double f, double g);
    virtual Intervald hgg(Intervald f, Intervald g);

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);
#endif

    virtual Intervald proc(const Box<double>& );
    virtual Box3d grad(const Box<double>& );
    virtual IMatrix3d hess(const Box<double>& );

	virtual void procq(const gmVector3 & x, double* q);
};

#endif

