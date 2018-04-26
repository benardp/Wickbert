/** 
 * @file Blend.h 
 * Definition of the abstract Blend class.
 * 
 * @author William Nagel (for CS 497).
 * @date Fall 2000
 */

#ifndef BLEND_H
#define BLEND_H

#include "BinaryOp.h"

/** Blend is an abstract class designed to support a variety of blending
 * operation.
 *
 * The main functionality of blend is the definition of the method h(f,g) which
 * is a real function of two real variables f and g which are the values of the
 * two implicits to blend.
 */
class Blend : public BinaryOp
{
  private:
    void init(Implicit*,Implicit*);

  public:
    double m_r1, m_r2;  // blend radii used only for local blends -jh

    /// Constructor
	Blend(Implicit *f = NULL, Implicit *g = NULL, double r1 = 1.0, double r2 = 1.0) :
		BinaryOp(f,g), m_r1(r1), m_r2(r2) { }

    /** The function h merges the value of the two implicits to blend.
     * New blends need only implement h and its derivatives.
     *
     * \parameter f,g The values of the two implicits.
     * \return The value of the blend of the two implicits.
     */
    virtual double h(double f, double g) = 0;
    virtual Intervald h(Intervald f, Intervald g) = 0;

    /** dh/df
     * Defaults to a numerical approximation.
     */
    virtual double hf(double f, double g) 
      { return (h(f+m_epsilon,g) - h(f,g))/m_epsilon; }
    virtual Intervald hf(Intervald f, Intervald g) 
      { return (h(f+Intervald(m_epsilon),g) - h(f,g))/Intervald(m_epsilon); }

    /** dh/dg
     * Defaults to a numerical approximation.
     */
    virtual double hg(double f, double g) 
      { return (h(f,g+m_epsilon) - h(f,g))/m_epsilon; }
    virtual Intervald hg(Intervald f, Intervald g) 
      { return (h(f,g+Intervald(m_epsilon)) - h(f,g))/Intervald(m_epsilon); }

    /** d^2h/df^2
     * Defaults to a numerical approximation.
     */
    virtual double hff(double f, double g) 
      { return (hf(f+m_epsilon,g) - hf(f,g))/m_epsilon; }
    virtual Intervald hff(Intervald f, Intervald g) 
      { return (hf(f+Intervald(m_epsilon),g) - hf(f,g))/Intervald(m_epsilon); }

    /** d^2h/dg^2
     * Defaults to a numerical approximation.
     */
    virtual double hgg(double f, double g) 
      { return (hf(f,g+m_epsilon) - hf(f,g))/m_epsilon; }
    virtual Intervald hgg(Intervald f, Intervald g) 
      { return (hf(f,g+Intervald(m_epsilon)) - hf(f,g))/Intervald(m_epsilon); }

    /** d^2h/dfdg
     * Defaults to a numerical approximation.
     */
    virtual double hfg(double f, double g) 
      { return (hg(f+m_epsilon,g) - hg(f,g))/m_epsilon; }
    virtual Intervald hfg(Intervald f, Intervald g) 
      { return (hg(f+Intervald(m_epsilon),g) - hg(f,g))/Intervald(m_epsilon); }

    /** dh/dr1 The derivative of h wrt the first radius.
     * Defaults to zero since we have no way of knowing.
     */
    virtual double hr1(double f, double g) { return 0.0; }
    virtual Intervald hr1(Intervald f, Intervald g) { return Intervald(0.0); }

    /** dh/dr2 The derivative of h wrt the second radius.
     * Defaults to zero since we have no way of knowing.
     */
    virtual double hr2(double f, double g) { return 0.0; }
    virtual Intervald hr2(Intervald f, Intervald g) { return Intervald(0.0); }

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);
#endif

    virtual Intervald proc(const Box<double>& );
    virtual Box3d grad(const Box<double>& );
    virtual IMatrix3d hess(const Box<double>& );

    virtual void _setq(double* q);
    virtual void getq(double* q);
    virtual void procq(const gmVector3 & x, double* q);
    virtual void getqname(char** name);
    virtual unsigned int qlen(void);
};

#endif 

