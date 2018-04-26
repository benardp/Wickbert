/**
 * Implementation of the RFunction class.
 *
 * @author Steve Zelinka
 */

#include "RFunction.h"

RFunction::RFunction(double cont, double sign)
{
	m_cont = cont;
	m_sign = sign;

	new SurfParamDouble(this, &m_cont, cont, "Continuity", "R-Function Continuity exponent", "R-Function Continuity exponent");
	new SurfParamDouble(this, &m_sign, sign, "Sign", "R-Function Sign", "R-Function Sign");
}

double RFunction::h(double f, double g)
{
  double fg = sqrt(f*f + g*g);
  double h = f + g + m_sign*fg;
  h *= pow(fg,m_cont);
  return h;
}

Intervald RFunction::h(Intervald f, Intervald g)
{
  Intervald fg = (f.squared() + g.squared()).sqrt();
  Intervald h  = f + g + Intervald(m_sign)*fg;
  h *= fg.pow(Intervald(m_cont));
  return h;
}

double RFunction::hf(double f, double g)
{
#if 0		//	should write the equation as a polynomial in r then take the derivative, to avoid dividing by zero in the intermediate steps		-Xinlai
  /*
   * h = (f + g + m_sign*(sqrt(f^2+g^2)))*(f^2+g^2)^(m_cont/2).
   * hf = (f^2+g^2)^(m_cont/2)*d(f + g + m_sign*(sqrt(f^2+g^2)))/df +
   *      (f + g + m_sign*(sqrt(f^2+g^2)))*d(f^2+g^2)^(m_cont/2)/df
   *    = (f^2+g^2)^(m_cont/2)*(1 + m_sign*f/sqrt(f^2+g^2) +
   *      (f + g + m_sign*(sqrt(f^2+g^2)))*m_cont*f*(f^2+g^2)^(m_cont/2 - 1).
   */

  double fg = sqrt(f*f + g*g);
  double hf = pow(fg,m_cont)*(1.0 + m_sign*f/fg) + 
             (f + g + m_sign*fg)*m_cont*f*pow(fg,m_cont - 2.0);
  return hf;
#endif
  /**
   *	let r = sqrt(f^2 + g^2)
   *	h = (f + g + m_sign * r) * r^m
   *	  = f * r^m + g * r^m + m_sign * r^{m+1}
   *	hf = r^m + m(f^2 + f*g)r^{m-2} + m_sign * f * r^{m-1}
   */
  double r = sqrt(f * f + g * g);
  double hf = pow(r, m_cont) + m_cont * (f * f + f * g) * pow(r, m_cont - 2) + m_sign * (m_cont + 1) * f * pow(r, m_cont - 1);
  return hf;
}

Intervald RFunction::hf(Intervald f, Intervald g)
{
#if 0	//	should write the equation as a polynomial in r then take the derivative, to avoid dividing by zero in the intermediate steps		-Xinlai
  /*
   * h = (f + g + m_sign*(sqrt(f^2+g^2)))*(f^2+g^2)^(m_cont/2).
   * hf = (f^2+g^2)^(m_cont/2)*d(f + g + m_sign*(sqrt(f^2+g^2)))/df +
   *      (f + g + m_sign*(sqrt(f^2+g^2)))*d(f^2+g^2)^(m_cont/2)/df
   *    = (f^2+g^2)^(m_cont/2)*(1 + m_sign*f/sqrt(f^2+g^2) +
   *      (f + g + m_sign*(sqrt(f^2+g^2)))*m_cont*f*(f^2+g^2)^(m_cont/2 - 1).
   */

  Intervald fg = (f.squared() + g.squared()).sqrt();
  Intervald hf = fg.pow(Intervald(m_cont)) * 
                 (Intervald(1.0) + Intervald(m_sign) * f / fg) +
                 (f + g + Intervald(m_sign)*fg) * Intervald(m_cont) * f * 
      fg.pow(Intervald(m_cont-2.0));
#endif
	/**
   *	let r = sqrt(f^2 + g^2)
   *	h = (f + g + m_sign * r) * r^m
   *	  = f * r^m + g * r^m + m_sign * r^{m+1}
   *	hf = r^m + m(f^2 + f*g)r^{m-2} + m_sign * f * r^{m-1}
   */
  Intervald r = (f.squared() + g.squared()).sqrt();
  Intervald hf = r.pow(m_cont) + Intervald(m_cont) * Intervald(f * f + f * g) * r.pow(m_cont - 2) + Intervald(m_sign * (m_cont + 1) * f) * r.pow(m_cont - 1);
  return hf;
}

double RFunction::hg(double f, double g)
{
#if 0	//	should write the equation as a polynomial in r then take the derivative, to avoid dividing by zero in the intermediate steps		-Xinlai
  /*
   * h = (f + g + m_sign*(sqrt(f^2+g^2)))*(f^2+g^2)^(m_cont/2).
   * hg = (f^2+g^2)^(m_cont/2)*d(f + g + m_sign*(sqrt(f^2+g^2)))/dg +
   *      (f + g + m_sign*(sqrt(f^2+g^2)))*d(f^2+g^2)^(m_cont/2)/dg
   *    = (f^2+g^2)^(m_cont/2)*(1 + m_sign*g/sqrt(f^2+g^2) +
   *      (f + g + m_sign*(sqrt(f^2+g^2)))*m_cont*g*(f^2+g^2)^(m_cont/2 - 1).
   */
  double fg = sqrt(f*f + g*g);
  double hg = pow(fg,m_cont)*(1.0 + m_sign*g/fg) + 
              (f + g + m_sign*fg)*m_cont*g*pow(fg,m_cont - 2.0);
#endif
	/**
   *	let r = sqrt(f^2 + g^2)
   *	h = (f + g + m_sign * r) * r^m
   *	  = f * r^m + g * r^m + m_sign * r^{m+1}
   *	hg = r^m + m(g^2 + f*g)r^{m-2} + m_sign * g * r^{m-1}
   */
  double r = sqrt(f * f + g * g);
  double hg = pow(r, m_cont) + m_cont * (g * g + g * f) * pow(r, m_cont - 2) + m_sign * (m_cont + 1) * g * pow(r, m_cont - 1);
  return hg;
}

Intervald RFunction::hg(Intervald f, Intervald g)
{
	#if 0	//	should write the equation as a polynomial in r then take the derivative, to avoid dividing by zero in the intermediate steps		-Xinlai

  /*
   * h = (f + g + m_sign*(sqrt(f^2+g^2)))*(f^2+g^2)^(m_cont/2).
   * hg = (f^2+g^2)^(m_cont/2)*d(f + g + m_sign*(sqrt(f^2+g^2)))/dg +
   *      (f + g + m_sign*(sqrt(f^2+g^2)))*d(f^2+g^2)^(m_cont/2)/dg
   *    = (f^2+g^2)^(m_cont/2)*(1 + m_sign*g/sqrt(f^2+g^2) +
   *      (f + g + m_sign*(sqrt(f^2+g^2)))*m_cont*g*(f^2+g^2)^(m_cont/2 - 1).
   */
  Intervald fg = (f.squared() + g.squared()).sqrt();
  Intervald hg = fg.pow(Intervald(m_cont)) * 
                 (Intervald(1.0) + Intervald(m_sign) * g / fg) +
                 (f + g + Intervald(m_sign)*fg) * Intervald(m_cont) * g *
      fg.pow(Intervald(m_cont-2.0));
#endif
  /**
   *	let r = sqrt(f^2 + g^2)
   *	h = (f + g + m_sign * r) * r^m
   *	  = f * r^m + g * r^m + m_sign * r^{m+1}
   *	hg = r^m + m(g^2 + f*g)r^{m-2} + m_sign * g * r^{m-1}
   */
  Intervald r = (f.squared() + g.squared()).sqrt();
  Intervald hg = r.pow(m_cont) + Intervald(m_cont) * Intervald(g * g + g * f) * r.pow(m_cont - 2) + Intervald(m_sign * (m_cont + 1) * g) * r.pow(m_cont - 1);
  return hg;
}

double RFunction::hff(double f, double g)
{
  /*
   * fg = sqrt(f^2 + g^2)
   * dfg/df = f/fg
   * hf = fg^m_cont*(1 + m_sign*f/fg) +
   *      (f + g + m_sign*fg)*m_cont*f*fg^(m_cont - 2)
   *    = fg^m_cont + m_sign*f*fg^(m_cont-1) +
   *      (f + g + m_sign*fg)*m_cont*f*fg^(m_cont - 2)
   * hff = m_cont*fg^(m_cont-1)*dfg/df +
   *       m_sign*f*(m_cont-1)*fg^(m_cont-2)*dfg/df +
   *       m_sign*fg^(m_cont-1) +
   *       (1 + m_sign*dfg/df)*m_cont*f*fg^(m_cont - 2) +
   *       (f + g + m_sign*fg)*m_cont*(f*(m_cont - 2)*fg^(m_cont - 3)*dfg/df + 
   *       fg^(m_cont - 2))
   *     = 2*m_cont*f*fg^(m_cont-2) +
   *       m_sign*(m_cont-1)*f^2*fg^(m_cont-3) +
   *       m_sign*fg^(m_cont-1) +
   *       m_cont*m_sign*f^2*fg^(m_cont - 3) +
   *       m_cont*(f + g + m_sign*fg)*((m_cont - 2)*f^2*fg^(m_cont - 4) + 
   *       fg^(m_cont - 2))
   */

  double fg = sqrt(f*f + g*g);
  double hff;
  
  hff = 2.0*m_cont*f*pow(fg,m_cont-2.0);
  hff += m_sign*(m_cont-1.0)*f*f*pow(fg,m_cont-3.0);
  hff += m_sign*pow(fg,m_cont-1.0);
  hff += m_cont*m_sign*f*f*pow(fg,m_cont-3.0);
  hff += m_cont*(f + g + m_sign*fg)*((m_cont - 2.0)*f*f*pow(fg,m_cont-4.0) + 
         pow(fg,m_cont-2.0));

  return hff;
}

Intervald RFunction::hff(Intervald f, Intervald g)
{
  /*
   * fg = sqrt(f^2 + g^2)
   * dfg/df = f/fg
   * hf = fg^m_cont*(1 + m_sign*f/fg) +
   *      (f + g + m_sign*fg)*m_cont*f*fg^(m_cont - 2)
   *    = fg^m_cont + m_sign*f*fg^(m_cont-1) +
   *      (f + g + m_sign*fg)*m_cont*f*fg^(m_cont - 2)
   * hff = m_cont*fg^(m_cont-1)*dfg/df +
   *       m_sign*f*(m_cont-1)*fg^(m_cont-2)*dfg/df +
   *       m_sign*fg^(m_cont-1) +
   *       (1 + m_sign*dfg/df)*m_cont*f*fg^(m_cont - 2) +
   *       (f + g + m_sign*fg)*m_cont*(f*(m_cont - 2)*fg^(m_cont - 3)*dfg/df + 
   *       fg^(m_cont - 2))
   *     = 2*m_cont*f*fg^(m_cont-2) +
   *       m_sign*(m_cont-1)*f^2*fg^(m_cont-3) +
   *       m_sign*fg^(m_cont-1) +
   *       m_cont*m_sign*f^2*fg^(m_cont - 3) +
   *       m_cont*(f + g + m_sign*fg)*((m_cont - 2)*f^2*fg^(m_cont - 4) + 
   *       fg^(m_cont - 2))
   */

  Intervald fg = (f.squared() + g.squared()).sqrt();
  Intervald hff;
  
  hff =  Intervald(2.0)*Intervald(m_cont)*f*fg.pow(Intervald(m_cont-2.0));
  hff += Intervald(m_sign)*Intervald(m_cont-1.0)*f.squared()*
         fg.pow(Intervald(m_cont-3.0));
  hff += Intervald(m_sign)*fg.pow(Intervald(m_cont-1.0));
  hff += Intervald(m_cont)*Intervald(m_sign)*f.squared()*
         fg.pow(Intervald(m_cont-3.0));
  hff += Intervald(m_cont)*(f + g + Intervald(m_sign)*fg)*
         (Intervald(m_cont-2.0)*f.squared()*fg.pow(Intervald(m_cont-4.0)) + 
   fg.pow(Intervald(m_cont-2.0)));

  return hff;
}

double RFunction::hgg(double f, double g)
{
  /*
   * fg = sqrt(f^2 + g^2)
   * dfg/dg = g/fg
   * hg = fg^m_cont*(1 + m_sign*g/fg) +
   *      (f + g + m_sign*fg)*m_cont*g*fg^(m_cont - 2)
   *    = fg^m_cont + m_sign*g*fg^(m_cont-1) +
   *      (f + g + m_sign*fg)*m_cont*g*fg^(m_cont - 2)
   * hgg = m_cont*fg^(m_cont-1)*dfg/dg +
   *       m_sign*g*(m_cont-1)*fg^(m_cont-2)*dfg/dg +
   *       m_sign*fg^(m_cont-1) +
   *       (1 + m_sign*dfg/dg)*m_cont*g*fg^(m_cont - 2) +
   *       (f + g + m_sign*fg)*m_cont*(g*(m_cont - 2)*fg^(m_cont - 3)*dfg/dg + 
   *       fg^(m_cont - 2))
   *     = 2*m_cont*g*fg^(m_cont-2) +
   *       m_sign*(m_cont-1)*g^2*fg^(m_cont-3) +
   *       m_sign*fg^(m_cont-1) +
   *       m_cont*m_sign*g^2*fg^(m_cont - 3) +
   *       m_cont*(f + g + m_sign*fg)*((m_cont - 2)*g^2*fg^(m_cont - 4) + 
   *       fg^(m_cont - 2))
   */
     
  double fg = sqrt(f*f + g*g);
  double hgg;
  
  hgg = 2.0*m_cont*g*pow(fg,m_cont-2.0);
  hgg += m_sign*(m_cont-1.0)*g*g*pow(fg,m_cont-3.0);
  hgg += m_sign*pow(fg,m_cont-1.0);
  hgg += m_cont*m_sign*g*g*pow(fg,m_cont-3.0);
  hgg += m_cont*(f + g + m_sign*fg)*((m_cont - 2.0)*g*g*pow(fg,m_cont-4.0) + 
         pow(fg,m_cont-2.0));

  return hgg;
}

Intervald RFunction::hgg(Intervald f, Intervald g)
{
  /*
   * fg = sqrt(f^2 + g^2)
   * dfg/dg = g/fg
   * hg = fg^m_cont*(1 + m_sign*g/fg) +
   *      (f + g + m_sign*fg)*m_cont*g*fg^(m_cont - 2)
   *    = fg^m_cont + m_sign*g*fg^(m_cont-1) +
   *      (f + g + m_sign*fg)*m_cont*g*fg^(m_cont - 2)
   * hgg = m_cont*fg^(m_cont-1)*dfg/dg +
   *       m_sign*g*(m_cont-1)*fg^(m_cont-2)*dfg/dg +
   *       m_sign*fg^(m_cont-1) +
   *       (1 + m_sign*dfg/dg)*m_cont*g*fg^(m_cont - 2) +
   *       (f + g + m_sign*fg)*m_cont*(g*(m_cont - 2)*fg^(m_cont - 3)*dfg/dg + 
   *       fg^(m_cont - 2))
   *     = 2*m_cont*g*fg^(m_cont-2) +
   *       m_sign*(m_cont-1)*g^2*fg^(m_cont-3) +
   *       m_sign*fg^(m_cont-1) +
   *       m_cont*m_sign*g^2*fg^(m_cont - 3) +
   *       m_cont*(f + g + m_sign*fg)*((m_cont - 2)*g^2*fg^(m_cont - 4) + 
   *       fg^(m_cont - 2))
   */
     
  Intervald fg = (f.squared() + g.squared()).sqrt();
  Intervald hgg;
  
  hgg =  Intervald(2.0)*Intervald(m_cont)*g*fg.pow(Intervald(m_cont-2.0));
  hgg += Intervald(m_sign)*(m_cont-1.0)*g.squared()*
         fg.pow(Intervald(m_cont-3.0));
  hgg += Intervald(m_sign)*fg.pow(Intervald(m_cont-1.0));
  hgg += Intervald(m_cont)*Intervald(m_sign)*g.squared()*
         fg.pow(Intervald(m_cont-3.0));
  hgg += Intervald(m_cont)*(f + g + Intervald(m_sign)*fg)*
         (Intervald(m_cont-2.0)*g.squared()*fg.pow(Intervald(m_cont-4.0)) + 
         fg.pow(Intervald(m_cont-2.0)));

  return hgg;
}


double RFunction::hfg(double f, double g)
{
//	old code derivation error
#if 0
  /*
   * fg = sqrt(f^2 + g^2)
   * dfg/dg = g/fg
   * hf = fg^m_cont*(1 + m_sign*f/fg) +
   *      (f + g + m_sign*fg)*m_cont*f*fg^(m_cont - 2)
   *    = fg^m_cont + m_sign*f*fg^(m_cont-1) +
   *      (f + g + m_sign*fg)*m_cont*f*fg^(m_cont - 2)
   * hfg = m_cont*fg^(m_cont-1)*dfg/dg +
   *       m_sign*f*(m_cont-1)*fg^(m_cont-2)*dfg/dg +
   *       m_sign*fg^(m_cont-1) +					//	derivation error here, it seems the original author is trying to take the derivative of m_sign*f*fg^(m_cont-2) wrt f
   *       (1 + m_sign*dfg/dg)*m_cont*f*fg^(m_cont - 2) +
   *       (f + g + m_sign*fg)*m_cont*(f*(m_cont - 2)*fg^(m_cont - 3)*dfg/dg + 
   *       fg^(m_cont - 2))
   *     = m_sign*fg^(m_cont-1.0) +
   *       m_cont*(f+g)*fg^(m_cont-2.0) +
   *       m_sign*(2.0*m_cont-1.0)*f*g*fg^(m_cont-3.0) +
   *       m_cont*(f + g + m_sign*fg)*((m_cont - 2)*f*g*fg^(m_cont - 4.0) + 
   *       pow(fg,m_cont - 2.0));
   *
   *     = hgf
   */

  double fg = sqrt(f*f + g*g);
  double hfg;
  hfg = m_cont * pow(fg, m_cont - 1.0);
  hfg += m_cont*(f+g)*pow(fg,m_cont-2.0);
  hfg += m_sign*(2.0*m_cont-1.0)*f*g*pow(fg,m_cont-3.0);
  hfg += m_cont*(f + g + m_sign*fg)*((m_cont - 2.0)*f*g*pow(fg,m_cont-4.0) + 
         pow(fg,m_cont-2.0));
#endif
  /*
   *	fg = sqrt(f^2 + g^2)
   *	dfg/dg = g/fg
   *	hf = fg^m_cont(1 + m_sign*f/fg) +
   *		 (f + g + m_sign * fg) * m_cont * f * fg^(m_cont - 2)
   *		= fg^m_cont + m_sign*f*fg^(m_cont - 1) +
   *		  (f + g + m_sign * fg) * m_cont * f * fg^(m_cont-2)
   *	hfg = m_cont * fg^(m_cont - 1) * dfg/dg +
   *		  m_sign*f*(m_cont - 1)*fg^(m_cont - 2) * dfg/dg +
   *		  (1+m_sign*dfg/dg) * m_cont * f * fg^(m_cont-2) +
   *		  (f + g + m_sign * fg) * m_cont * f * (m_cont-2) * fg^(m_cont-3) * dfg/dg
   *		= m_cont * g * fg^(m_cont - 2) +
   *		  m_sign * f * g * (m_cont - 1) * fg^(m_cont - 3) +
   *		  m_cont * f * fg^(m_cont - 2) + m_sign * m_cont * f * g * fg^(m_cont - 3) +
   *		  (f + g + m_sign * fg) * m_cont * f * g * (m_cont - 2) * fg^(m_cont - 4)
   *		= m_cont * (g + f) * fg^(m_cont - 2) +
   *		  m_sign * f * g * (2*m_cont-1)*fg^(m_cont-3) +
   *		  (f + g + m_sign * fg) * m_cont * (m_cont - 2) * f * g * fg^(m_cont-4)
   */
  double fg = sqrt(f*f + g*g);
  double hfg;
  hfg = m_cont*(f+g)*pow(fg,m_cont-2.0);
  hfg += m_sign*(2.0*m_cont-1.0)*f*g*pow(fg,m_cont-3.0);
  hfg += m_cont*(f + g + m_sign*fg)*(m_cont - 2.0)*f*g*pow(fg,m_cont-4.0);
  return hfg;
}

Intervald RFunction::hfg(Intervald f, Intervald g)
{
//	old code has derivation error
#if 0
  /*
   * fg = sqrt(f^2 + g^2)
   * dfg/dg = g/fg
   * hf = fg^m_cont*(1 + m_sign*f/fg) +
   *      (f + g + m_sign*fg)*m_cont*f*fg^(m_cont - 2)
   *    = fg^m_cont + m_sign*f*fg^(m_cont-1) +
   *      (f + g + m_sign*fg)*m_cont*f*fg^(m_cont - 2)
   * hfg = m_cont*fg^(m_cont-1)*dfg/dg +
   *       m_sign*f*(m_cont-1)*fg^(m_cont-2)*dfg/dg +
   *       m_sign*fg^(m_cont-1) +
   *       (1 + m_sign*dfg/dg)*m_cont*f*fg^(m_cont - 2) +
   *       (f + g + m_sign*fg)*m_cont*(f*(m_cont - 2)*fg^(m_cont - 3)*dfg/dg + 
   *       fg^(m_cont - 2))
   *     = m_sign*fg^(m_cont-1.0) +
   *       m_cont*(f+g)*fg^(m_cont-2.0) +
   *       m_sign*(2.0*m_cont-1.0)*f*g*fg^(m_cont-3.0) +
   *       m_cont*(f + g + m_sign*fg)*((m_cont - 2)*f*g*fg^(m_cont - 4.0) + 
   *       pow(fg,m_cont - 2.0));
   *
   *     = hgf
   */

  Intervald fg = (f.squared() + g.squared()).sqrt();
  Intervald hfg;
  
  hfg =  Intervald(m_sign)*fg.pow(Intervald(m_cont-1.0));
  hfg += Intervald(m_cont)*(f+g)*fg.pow(Intervald(m_cont-2.0));
  hfg += Intervald(m_sign)*Intervald(2.0*m_cont-1.0)*f*g*
         fg.pow(Intervald(m_cont-3.0));
  hfg += Intervald(m_cont)*(f + g + Intervald(m_sign)*fg)*
         (Intervald(m_cont-2.0)*f*g*fg.pow(Intervald(m_cont-4.0)) + 
         fg.pow(Intervald(m_cont-2.0)));
#endif
/*
   *	fg = sqrt(f^2 + g^2)
   *	dfg/dg = g/fg
   *	hf = fg^m_cont(1 + m_sign*f/fg) +
   *		 (f + g + m_sign * fg) * m_cont * f * fg^(m_cont - 2)
   *		= fg^m_cont + m_sign*f*fg^(m_cont - 1) +
   *		  (f + g + m_sign * fg) * m_cont * f * fg^(m_cont-2)
   *	hfg = m_cont * fg^(m_cont - 1) * dfg/dg +
   *		  m_sign*f*(m_cont - 1)*fg^(m_cont - 2) * dfg/dg +
   *		  (1+m_sign*dfg/dg) * m_cont * f * fg^(m_cont-2) +
   *		  (f + g + m_sign * fg) * m_cont * f * (m_cont-2) * fg^(m_cont-3) * dfg/dg
   *		= m_cont * g * fg^(m_cont - 2) +
   *		  m_sign * f * g * (m_cont - 1) * fg^(m_cont - 3) +
   *		  m_cont * f * fg^(m_cont - 2) + m_sign * m_cont * f * g * fg^(m_cont - 3) +
   *		  (f + g + m_sign * fg) * m_cont * f * g * (m_cont - 2) * fg^(m_cont - 4)
   *		= m_cont * (g + f) * fg^(m_cont - 2) +
   *		  m_sign * f * g * (2*m_cont-1)*fg^(m_cont-3) +
   *		  (f + g + m_sign * fg) * m_cont * (m_cont - 2) * f * g * fg^(m_cont-4)
   */
  Intervald fg = (f.squared() + g.squared()).sqrt();
  Intervald hfg;
  
  hfg = Intervald(m_cont)*(f+g)*fg.pow(Intervald(m_cont-2.0));
  hfg += Intervald(m_sign)*Intervald(2.0*m_cont-1.0)*f*g*
         fg.pow(Intervald(m_cont-3.0));
  hfg += Intervald(m_cont)*(f + g + Intervald(m_sign)*fg)*
         Intervald(m_cont-2.0)*f*g*fg.pow(Intervald(m_cont-4.0));
  return hfg;
}


#ifndef INTERVAL_EVAL_ONLY
double RFunction::proc(const gmVector3 & x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    return h(m_f->proc(x),m_g->proc(x));
  else
    return 0.0;
}

gmVector3 RFunction::grad(const gmVector3 & x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    {
      double fx = m_f->proc(x); 
      double gx = m_g->proc(x);
      return hf(fx,gx) * m_f->grad(x) + hg(fx,gx) * m_g->grad(x);
    }
  else
    return gmVector3();
}

gmMatrix3 RFunction::hess(const gmVector3 & x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    {
      double fx = m_f->proc(x); 
      double gx = m_g->proc(x);
      gmVector3 dfx = m_f->grad(x);
      gmVector3 dgx = m_g->grad(x);

      return hff(fx,gx) * outer(dfx,dfx) + 
             hfg(fx,gx) * outer(dfx,dgx) +
             hfg(fx,gx) * outer(dgx,dfx) + 
             hgg(fx,gx) * outer(dgx,dgx) +
             hf(fx,gx)  * m_f->hess(x)   + 
             hg(fx,gx)  * m_g->hess(x);
    }
  else
    return gmMatrix3();
}
#endif

Intervald RFunction::proc(const Box<double>& x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    return h(m_f->proc(x),m_g->proc(x));
  else
    return Intervald(0.0);
}

Box3d RFunction::grad(const Box<double>& x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    {
      Intervald fx = m_f->proc(x);
      Intervald gx = m_g->proc(x);
      return hf(fx,gx) * m_f->grad(x) + hg(fx,gx) * m_g->grad(x);
    }
  else
    return Box3d(0.0);
}

IMatrix3d RFunction::hess(const Box<double>& x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    {
      Intervald fx = m_f->proc(x); 
      Intervald gx = m_g->proc(x);
      Box3d dfx = m_f->grad(x);
      Box3d dgx = m_g->grad(x);

      return hff(fx,gx) * outer(dfx,dfx) + 
             hfg(fx,gx) * outer(dfx,dgx) +
             hfg(fx,gx) * outer(dgx,dfx) + 
             hgg(fx,gx) * outer(dgx,dgx) +
             hf(fx,gx)  * m_f->hess(x)   + 
             hg(fx,gx)  * m_g->hess(x);
    }
  else
    return IMatrix3d(0.0);
}

void RFunction::procq(const gmVector3 & x, double* q)
{
	if( m_f != NULL && m_g != NULL)
	{
		m_f->procq(x, q);
		m_g->procq(x, &q[m_f->qlen()]);
		double fx = m_f->proc(x);
		double gx = m_g->proc(x);
		double dhdf = hf(fx, gx);
		unsigned int i = 0;
		for( i = 0; i < m_f->qlen(); i++)
			q[i] *= dhdf;
		double dhdg = hg(fx, gx);
		for(; i < m_f->qlen() + m_g->qlen(); i++)
			q[i] *= dhdg;
	}

}

