/**
 * @file Blinn.h 
 */

#ifndef __BLINN_H
#define __BLINN_H

#include "Blob.h"
#include "Surface/Implicit/Algebraic/Quadric.h"

/**
 * Class Blinn -- Blob (Gaussian) with
 * an arbitrary quadric radius function
 */
class Blinn : public Blob
{
public:
	double m_b;

    void init(Quadric *f, double b, double r);

    /// Class Constructors                         
    Blinn();
    Blinn(Quadric *f);
    Blinn(Quadric *f, double b);
    Blinn(Quadric *f, double b, double r);
    
    /// Kernel functions
    virtual double kernel(double r2);
    virtual double dkernel(double r2);
    virtual double d2kernel(double r2);

    /// Interval versions of Kernel functions
    virtual Intervald kernel(Intervald r2);
    virtual Intervald dkernel(Intervald r2);
    virtual Intervald d2kernel(Intervald r2);

    /// Accessor functions                      
    virtual unsigned int  qlen(void);
    virtual void _setq(double *q);
    virtual void getq(double *q);
    virtual void getqname(char **qn);

    virtual void procq(const gmVector3& v, double* dfdq);

    MAKE_NAME();
};

#endif

