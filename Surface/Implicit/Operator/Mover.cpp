/**
 * @file Mover.cpp
 */
#include "Mover.h"
#include <stdlib.h>

// TNT includes
#include "tnt/tnt.h"
#include "tnt/vec.h"
#include "tnt/cholesky.h"
#include "tnt/lu.h"
#include "tnt/cmat.h"
#include "tnt/transv.h"
#include "tnt/trisolve.h"
#include "Surface/RandomStream.h"

REGISTER_IMPLICIT(Mover,"UnaryOp:Mover");

/**
 * Default constructor.
 */
Mover::Mover()
{
  m_f = NULL;
  m_o = gmVector3();
}

/**
 * Convenience function for the calculation of a LRP over a time interval
 * for a particular m_x.  We use the "current" and the "last" coefficients
 * to calculate an interval of the coefficients over a given time interval.
 * This is used for 4D critical point finding.
 */
Box3d Mover::m_olrp(Intervald t)
{
  double* m_oold = (double *) malloc(qlen() * sizeof(double));
  getqold(m_oold);
  Box3d retval;

  // m_o resides in q[0],q[1],q[2] = m_oold[0..2]
  for (int i = 0; i < 3; i++)
    retval[i] = t * (m_o[i] - m_oold[i]) + Intervald(m_oold[i]);

  delete m_oold;
  return retval;
}

gmVector3 Mover::dinv(const gmVector3 & x)
{
  return x - m_o;
}

Box3d Mover::dinv(Box3d x)
{
  return x - Box3d(m_o);
}

Box4d Mover::dinv(Box4d x)
{
  Box3d x1 = Box3d(x[0],x[1],x[2]);
  Box3d x2 = x1 - m_olrp(x[3]);
  return Box4d(x2[0],x2[1],x2[2],x[3]);
}

gmMatrix3 Mover::dinvjac(const gmVector3 & x)
{
  return gmMatrix3::identity();
}

IMatrix3d Mover::dinvjac(const Box<double>& x)
{
  return IMatrix3d::identity();
}

#ifndef INTERVAL_EVAL_ONLY
double Mover::proc(const gmVector3 & x)
{
  return (m_f ? m_f->proc(dinv(x)) : 0.0);
}

/**
 * grad of f(dinv(x)) = grad f(dinv(x)) * dinvjac(x).
 */
gmVector3 Mover::grad(const gmVector3 & x)
{
  if (m_f)
    {
      gmMatrix3 j = dinvjac(x);
      return j.transpose() * m_f->grad(dinv(x));
    }
  else
    return gmVector3();
}

/** Need to fix hess (does nothing now)
 */
gmMatrix3 Mover::hess(const gmVector3 & x)
{
  return (m_f ? m_f->hess(x) : gmMatrix3::identity());
}
#endif

Intervald Mover::proc(const Box<double>& x)
{
  Intervald retval(0.0);
  if (m_f)
    {
      if (x.size() == 3)
        retval = m_f->proc(dinv(Box3d(x)));
      else
        retval = m_f->proc(dinv(Box4d(x)));
    }
  return retval;
}

Box3d Mover::grad(const Box<double>& x)
{
  if (m_f)
    {
      IMatrix3d j = dinvjac(x);
      IMatrix3d jtrans = TNT::transpose(j);
      if (x.size() == 3)
        return jtrans * m_f->grad(dinv(Box3d(x)));
      else
        return jtrans * m_f->gradt(dinv(Box4d(x)));
    }
  else
    return Box3d(0.0);
}

IMatrix3d Mover::hess(const Box<double>& x)
{
  return (m_f ? m_f->hess(x) : IMatrix3d::identity());
}

/**
 * @todo Fix me, but may not be calculatable.
 */
Intervald Mover::proct(const Box<double>& x)
{
  return Intervald(0.0);
}

/**
 * @todo Fix me, but may not be calculatable.
 */
Box3d Mover::gradt(const Box<double>& x)
{
  return Box3d(0.0);
}

void Mover::procq(const gmVector3 & x, double *q)
{
  if (!m_f)
    return;

  /* df(x-o)/dq = df(x-o)/d(x-o) * d(x-o)/do
   */
  gmVector3 dfdo =-grad(x);

  q[0] = dfdo[0];
  q[1] = dfdo[1];
  q[2] = dfdo[2];
}

void Mover::getq(double *q)
{
  q[0] = m_o[0];
  q[1] = m_o[1];
  q[2] = m_o[2];
}

void Mover::_setq(double *q)
{
  if (m_o != gmVector3(q[0],q[1],q[2]))
  {
    m_o = gmVector3(q[0],q[1],q[2]);
    impdiff();
    m_o = gmVector3(0.0,0.0,0.0);
  }
}

unsigned int Mover::qlen() { return 3; }
unsigned int Mover::plen() { return 3; }

void Mover::getqname(char** qn)
{
  qn[0] = "X";
  qn[1] = "Y";
  qn[2] = "Z";
}

bool Mover::solve()
{
  if (!m_f) 
    return false;

  int i, j, iter, n = 10;
  double t,tnext;
  std::vector<gmVector3> p(n);  // vector of particle positions
  std::vector<gmVector3> pnext(n);  // vector of particle positions
  std::vector<gmVector3> p1(n);  // vector of goal particle positions
  std::vector<gmVector3> p0(n);  // vector of start particle positions
  TNT::Vector<double> b(n),x(n),f(n);
  TNT::Matrix<double> A(n,n);
  TNT::Vector<double> fqi(m_f->qlen()),fqj(m_f->qlen());
  TNT::Vector<double> qdot(m_f->qlen(),0.0);

  /* Place n random particles
   */
  RandomStream rs;
  for (i = 0; i < n; i++) 
    {
      p1[i] = gmVector3(rs.next(),rs.next(),rs.next());
      p0[i] = dinv(p1[i]);
      f[i] = m_f->proc(p[i]);
    }

  for (iter = 0; iter < 50; iter++) 
    {
      t = 0.02*(double)iter;
      tnext = 0.02*(double)(iter+1);

      /* Compute F_x^i . (D(pi) - pi)
       */
      for (i = 0; i < n; i++) 
        {
          p[i] = t*p1[i] + (1.0-t)*p0[i];
          pnext[i] = tnext*p1[i] + (1.0-tnext)*p0[i];

          b[i] = dot(m_f->grad(p[i]),pnext[i] - p[i]);

          /* Add feedback term.
          */
          // b[i] += phi*(m_f->proc(p[i]) - f[i]);
        }

      /* Compute Fqi.Fqj
       */
      for (j = 0; j < n; j++) 
        {
          m_f->procq(dinv(p[j]),fqj);
          for (i = 0; i < n; i++) 
            {
              m_f->procq(dinv(p[i]),fqi);
              A[j][i] = 0.0;
              for (unsigned k = 0; k < m_f->qlen(); k++)
                A[j][i] += fqi[k]*fqj[k];
            }
        }

      /* Solve.
       */
      TNT::Matrix<double> L(n,n);

      if (TNT::Cholesky_upper_factorization(A, L) != 0)
        return false;

      TNT::Vector<double> y = TNT::Lower_triangular_solve(L, b);
      TNT::Transpose_View<TNT::Matrix<double> > T = TNT::Transpose_view(L);
      x = TNT::Upper_triangular_solve(T, y);

      for (i = 0; i < n; i++) 
        {
          m_f->procq(dinv(p[i]),fqi);
          for (unsigned k = 0; k < m_f->qlen(); k++)
            qdot[k] -= x[i]*fqi[k];
        }

      /* Set new parameters in m_f.
       */
      TNT::Vector<double> fq(m_f->qlen());
      m_f->getq(fq);
      for (unsigned i = 0; i < m_f->qlen(); i++)
        qdot[i] *= 0.5;
      fq = fq + qdot;
      m_f->setq(fq);
    }

  return true;
}

/** Use implicit differentiation to propogate change in deformation
 * parameters with change in deformed object parameters.
 *
 * Deformed object given by F(D^-1(x;k);q) = F(x;k;q)
 * where k is the vector of deformation parameters.
 *
 * Need to compute deltaq = (dq/dk) deltak
 *
 * dqi/dkj = -Fkj/Fqi
 */
bool Mover::impdiff()
{
  int i;
  int n = m_f->qlen();      // # of variables to solve
  int nk = qlen();        // # of adapter parameters
  TNT::Vector<gmVector3> x(n);  // random points on which to solve
  TNT::Matrix<double> Fk(n,nk);  // Matrix of dF(xi)/dk
  TNT::Vector<double> k(nk);    // Deformation parameters
  double *kf = new double[nk];    // Temp store for def. params
  TNT::Vector<double> Fkk(n);    // Fk * k
  TNT::Matrix<double> Fq(n,n);    // Matrix of dF(xi)/dq
  double *qf = new double[n];    // Operand parameters


  /* Fill x with random points in [-1,1]^3
   */
  RandomStream rs;
  for (i = 0; i < n; i++)
    x[i] = gmVector3(rs.next(-1.0,1.0),rs.next(-1.0,1.0),rs.next(-1.0,1.0));

  /* Compute Fki = dF(x[i];k;q)/dk
   */
  for (i = 0; i < n; i++)
    procq(x[i],Fk[i]);

  /* Fill k with adapter parameters
   */
  getq(kf);
  for (i = 0; i < nk; i++)
    k[i] = kf[i];

  /* Compute Fqi = dF(x[i];k;q)/dq
   */
  for (i = 0; i < n; i++)
    m_f->procq(x[i],Fq[i]);

  Fkk = Fk*k;

#define DEBUG
#ifdef DEBUG
  std::cout << "Fk:" << Fk;
  std::cout << "k:" << k;
  std::cout << "Fq:" << Fq;
#endif

  //we should free the memory...
  delete [] kf;
  delete [] qf;


  /* Need to solve Fq*q = Fkk
   */
  TNT::Vector<int> ipiv;      // temporary used to store pivot info

  /* Factor Fk in place
   */
  if (TNT::LU_factor(Fq,ipiv))
    return false;

#ifdef DEBUG
  std::cout << "ipiv:" << ipiv;
#endif

  /* Solve for new operand parameters q
   */
  if (TNT::LU_solve(Fq,ipiv,Fkk))
    return false;

  m_f->getq(qf);
  for (i = 0; i < n; i++)
    qf[i] += Fkk[i];
  m_f->_setq(qf);

#ifdef DEBUG
  std::cout << "delta q:" << Fkk;
#endif

  return true;
}

