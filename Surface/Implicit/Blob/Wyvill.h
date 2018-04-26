/**
 * @file Wyvill.h
 */

#ifndef __WYVILL_H
#define __WYVILL_H

#include "Blob.h"
#include "Surface/Implicit/Algebraic/Quadric.h"

/**
 * Class Wyvill -- Polynomial Blob with
 * an arbitrary quadric radius function
 */
class Wyvill : public Blob
{ 
  private:
    gmVector3 grad1(const gmVector3 & v);
    gmVector3 grad2(const gmVector3 & v);
    gmVector3 grad3(const gmVector3 & v);
    void init(Quadric*,double); 

  public:
    /// Class Constructors
    Wyvill();
    Wyvill(Quadric *f); 
    Wyvill(Quadric *f, double r); 

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

