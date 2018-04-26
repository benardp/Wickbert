/** \file Algebraic.h An algebraic implicit surface class.
 * \author Jeff Decker
 * \date Fall 2000.
 */

#ifndef __ALGEBRAIC_H
#define __ALGEBRAIC_H

#include "Surface/Implicit/Implicit.h"

/** An algebraic implicit surface class.
 * This class implements arbitrary polynomials of three variables.
 */

class Algebraic : public Implicit
{
  private:
    void init(int);

  public:
    Algebraic();
    Algebraic(int degree);
    ~Algebraic();

	//overridden loading function to increase degree to match number of parameters
	virtual bool readImplicit(std::ifstream &file,bool verbose);

#ifndef INTERVAL_EVAL_ONLY    
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & v);
    virtual gmMatrix3 hess(const gmVector3 & v);
#endif

    virtual Intervald proc(const Box<double>& );
    virtual Box3d grad(const Box<double>& );
    virtual IMatrix3d hess(const Box<double>& );

    /// Degree of polynomial.
    int degree();
    void degree(int);

    /// Number of coefficients.
    int getNumCoef();

    static int coefficients(int);

    /// Return the ordinal of the x^i y^j z^k coefficient.
    int getCoefIndex(int i, int j, int k);

	void pretty(std::string &s);

    double getCoef(int i, int j, int k);
    void setCoef(int i, int j, int k, double a);

    virtual void procq(const gmVector3&, double*);
    virtual void procq(const Box<double>& , Intervald*);
    virtual void gradq(const Box<double>& , Intervald*, Intervald*, Intervald*);

    virtual void getq(double*);
    virtual void _setq(double*);
    virtual unsigned int qlen();
    virtual void getqname(char** qn);

    MAKE_NAME();

  private:
    /** Array of coefficients.
     * Should be accessed through getCoef and setCoef.
     */
    double *m_a;
    
    /// Linear interpolate - 4D critical points
    Intervald m_alrp(int,Intervald);

  protected:
    /// Degree.
    int m_d;

    /** Number of terms.
     * Initialized by calcNumCoef.
     * Accessed via getNumCoef.
     */
    int m_numCoef;

    /// Array storing the exponents of x in coefficient order.
    int *m_x;

    /// Array storing the exponents of y in coefficient order.
    int *m_y;

    /// Array storing the exponents of z in coefficient order.
    int *m_z;

    /// Cached array of powers of input value of x in coefficient order.
    double *m_xPow;

    /// Cached array of powers of input value of y in coefficient order.
    double *m_yPow;

    /// Cached array of powers of input value of z in coefficient order.
    double *m_zPow;

    void initPowerArrays();
    void initXYZPowerArrays(const gmVector3 & v);
    void initDerivCoefs(void);
    void calcNumCoef();

    double dx(const gmVector3 & v);
    double dy(const gmVector3 & v);
    double dz(const gmVector3 & v);

    double dx2(const gmVector3 & v); 
    double dy2(const gmVector3 & v);
    double dz2(const gmVector3 & v);

    double dxdy(const gmVector3 & v);
    double dxdz(const gmVector3 & v);
    double dydz(const gmVector3 & v);
};

#endif

