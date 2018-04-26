/** 
 * @file Blob.h A blob abstract class
 * @author John C. Hart
 * @date Fall 2000
 */

#ifndef __BLOB_H
#define __BLOB_H

#include "Surface/Implicit/Operator/UnaryOp.h"

/**
 * A pure virtual Blob class collecting methods common to all blobby models.
 *
 * Blob is a unary operator that assumes a quadric operand that returns
 * radius squared. A separate real kernel function of one variable is
 * defined by its instances that accepts the squared radius and returns
 * a quasi-Gaussian bump.
 */

class Blob : public UnaryOp
{
public:
    double m_r;  ///< The radius of blob object

    /// The kernel function defined by instances.
    virtual double kernel(double r2) = 0;

    /** The derivatives (wrt r2, not r !!!) of the kernel
     * \todo These should default to central differences.
     */
    virtual double dkernel(double r2) = 0;
    virtual double d2kernel(double r2) = 0;

    virtual Intervald kernel(Intervald r2) = 0;
    virtual Intervald dkernel(Intervald r2) = 0;
    virtual Intervald d2kernel(Intervald r2) = 0;

	virtual double proc(const gmVector3 & x) { return m_f ? kernel(m_f->proc(x)) : 0.0; }
	virtual gmVector3 grad(const gmVector3 & x) { return m_f ? dkernel(m_f->proc(x)) * m_f->grad(x) : gmVector3(); }

    /** Blob Hessian
     * Uses chain rule.
     * d2/dx2 blob(x) = d/dx dkernel(m_f(x))dm_f(x)/dx
     * = d2kernel(m_f(x))/dx2 (dm_f(x)/dx)^2 + dkernel(m_f(x))/dx d2m_f(x)/dx2
     */
	virtual gmMatrix3 hess(const gmVector3 & x) {
		return m_f ? d2kernel(m_f->proc(x)) * outer(m_f->grad(x),m_f->grad(x)) +
					 dkernel(m_f->proc(x)) * m_f->hess(x) : gmMatrix3(); }

	virtual Intervald proc(const Box<double>& x) { return m_f ? kernel(m_f->proc(x)) : Intervald(0.0); }
	virtual Box3d grad(const Box<double>& x)
    {
      Box3d result = Box3d(0.0);
      if (m_f != NULL) {
          Intervald dkern = dkernel(m_f->proc(x));
          result = m_f->grad(x);
          for (int i = 0; i < 3; i++)
            result[i] *= dkern;
      }
      return result;
    }

    /// Interval version of hess()
    virtual IMatrix3d hess(const Box<double>& x)
    {
      IMatrix3d result = IMatrix3d(0.0);
      if (m_f != NULL)
        {
          Intervald dkern = dkernel(m_f->proc(x));
          Intervald d2kern = d2kernel(m_f->proc(x));
          IMatrix3d hessres = m_f->hess(x);
          IMatrix3d outerres = outer(m_f->grad(x),m_f->grad(x));
          for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
              result[i][j] = (dkern * hessres[i][j]) + 
                             (d2kern * outerres[i][j]);
        }
      return result;
    }

    virtual const char ** getPixmapXPM(const int& size) const
    {
      if (size <= 16)
        return (const char **)blob_pixmap16;
      else if (size <= 32)
        return (const char **)blob_pixmap32;
      else
        return (const char **)blob_pixmap48;
    }
};

#endif

