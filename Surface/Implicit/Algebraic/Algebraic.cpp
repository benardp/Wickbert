/** \file Algebraic.cpp 
 * Implementation of the algebraic implicit surface class methods.
 * \author Jeff Decker (plus extensive comments by John C. Hart)
 * \date Fall 2000.
 */

#include "Algebraic.h"

/** Add Algebraic to the ImplicitFactory
 */
REGISTER_IMPLICIT(Algebraic,"Algebraic");

/**
 * Initialization routine called by constructors.  This function allows all
 * of the main object init code to be in one place.
 */
void Algebraic::init(int deg)
{
  m_d = deg;

  calcNumCoef();  // Should return 1
  m_a = new double[m_numCoef];
  m_x = new int[m_numCoef];
  m_y = new int[m_numCoef];
  m_z = new int[m_numCoef];
  m_xPow = new double[m_d+1];
  m_yPow = new double[m_d+1];
  m_zPow = new double[m_d+1];

  initPowerArrays();
  for (int i=0; i<m_numCoef; i++)
    m_a[i] = 0.0;
}

/** Constucts a constant Algebraic.
 */
Algebraic::Algebraic() : Implicit()
{
  init(0);
}

Algebraic::Algebraic(int deg) : Implicit()
{
  init(deg);
}

void Algebraic::degree(int deg)
{
  int i;

  double* m_a_back;
  m_a_back = m_a;
  int* m_x_back = m_x;
  int* m_y_back = m_y;
  int* m_z_back = m_z;
  double* m_xPow_back = m_xPow;
  double* m_yPow_back = m_yPow;
  double* m_zPow_back = m_zPow;

  int old_numCoef = m_numCoef;
  int old_m_d = m_d;

  m_d = deg;

  calcNumCoef();

  m_a = new double[m_numCoef];
  m_x = new int[m_numCoef];
  m_y = new int[m_numCoef];
  m_z = new int[m_numCoef];
  m_xPow = new double[m_d+1];
  m_yPow = new double[m_d+1];
  m_zPow = new double[m_d+1];

  for (i = 0; i < gmMin(old_numCoef, m_numCoef); i++)
    {
      m_a[i] = m_a_back[i];
      m_x[i] = m_x_back[i];
      m_y[i] = m_y_back[i];
      m_z[i] = m_z_back[i];
    }

  for (i = 0; i < gmMin(old_m_d+1, m_d+1); i++)
    {
      m_xPow[i] = m_xPow_back[i];
      m_yPow[i] = m_yPow_back[i];
      m_zPow[i] = m_zPow_back[i];
    }

  for (i = (int)gmMin(old_numCoef, m_numCoef); i < m_numCoef; i++)
    m_a[i] = 0.0;

  initPowerArrays();

  delete[] m_a_back;
  delete[] m_x_back;
  delete[] m_y_back;
  delete[] m_z_back;
  delete[] m_xPow_back;
  delete[] m_yPow_back;
  delete[] m_zPow_back;
}

Algebraic::~Algebraic()
{
  delete[] m_a;
  delete[] m_x;
  delete[] m_y;
  delete[] m_z;

  delete[] m_xPow;
  delete[] m_yPow;
  delete[] m_zPow;
}

//When loading, increase degree if number of parameters exceed qlen
bool Algebraic::readImplicit(std::ifstream &file,bool verbose)
{
	std::vector<double> params;
	readParameters(file,params);

	if(params.size() != qlen())
	{
		//params doesn't match qlen, change degree
		int newdegree = 0;
		while (Algebraic::coefficients(newdegree) < (int)params.size())
			newdegree++;

		degree(newdegree);
	}
	setq(params);

	return true;
}


/**
 * Convenience function for the calculation of a LRP over a time interval
 * for a particular m_a.  We use the "current" and the "last" coefficients
 * to calculate an interval of the coefficients over a given time interval.
 * This is used for 4D critical point finding.
 */
Intervald Algebraic::m_alrp(int i, Intervald t)
{
  double* m_aold = (double *) malloc(m_numCoef * sizeof(double));
  getqold(m_aold);
  Intervald retval = t * (m_a[i] - m_aold[i]) + Intervald(m_aold[i]);
  delete m_aold;
  return retval;
}

#ifndef INTERVAL_EVAL_ONLY
double Algebraic::proc(const gmVector3 & v)
{
  int i;
  double val = 0.0;

  initXYZPowerArrays(v);
  for (i=0; i<m_numCoef; i++)
    val += m_a[i] * m_xPow[m_x[i]] * m_yPow[m_y[i]] * m_zPow[m_z[i]];
  return val;
}

gmVector3 Algebraic::grad(const gmVector3 & v)
{
  return gmVector3(dx(v), dy(v), dz(v));
}

gmMatrix3 Algebraic::hess(const gmVector3 & v)
{
    double dfdxy = dxdy(v), dfdxz = dxdz(v), dfdyz = dydz(v);
    double dfdxx = dx2(v),  dfdyy = dy2(v),  dfdzz = dz2(v);

    gmMatrix3 hess(dfdxx,   dfdxy,  dfdxz,
                   dfdxy,   dfdyy,  dfdyz,
                   dfdxz,   dfdyz,  dfdzz);
    return hess;
}
#endif


/** proc of an interval-vector (aka a Box) X.
 * If X is a 4D box then
 * we assume q is lerp between qold and q.
 *
 * x[3] = [told,t]
 * m_alrp (i,x[3]) returns [(a_i)old,a_i]
 * 
 */
Intervald Algebraic::proc(const Box<double>& x)
{
  int i;
  Intervald val(0.0);

  if (x.size() == 3) {
    for (i=0; i<m_numCoef; i++) 
      val += m_a[i]*x[0].pow(m_x[i])*x[1].pow(m_y[i])*x[2].pow(m_z[i]);
  } else {
    for (i=0; i<m_numCoef; i++) 
      val += m_alrp(i,x[3])*x[0].pow(m_x[i])*x[1].pow(m_y[i])*x[2].pow(m_z[i]);
  }
  return val;
}

Box3d Algebraic::grad(const Box<double>& x)
{
  Box3d dx(0.0);
  Intervald m_atemp;

  for (int i = 0; i < m_numCoef; i++) 
    {
      m_atemp = ((x.size() == 3) ? Intervald(m_a[i]) : m_alrp(i,x[3]));
      if (m_x[i] >= 1) 
        dx[0] += m_atemp * Intervald(m_x[i]) *
                 (x[0].pow(m_x[i]-1)*x[1].pow(m_y[i])  *x[2].pow(m_z[i]));
      if (m_y[i] >= 1) 
        dx[1] += m_atemp * Intervald(m_y[i]) *
                 (x[0].pow(m_x[i])  *x[1].pow(m_y[i]-1)*x[2].pow(m_z[i]));
      if (m_z[i] >= 1) 
        dx[2] += m_atemp * Intervald(m_z[i]) *
                 (x[0].pow(m_x[i])  *x[1].pow(m_y[i])  *x[2].pow(m_z[i]-1));
    }
  return dx;
}

IMatrix3d Algebraic::hess(const Box<double>& x)
{
  IMatrix3d Vx;
  int i;

  // Calc lower triangle
  for (i = 0; i < m_numCoef; i++) 
    {
      if (m_x[i] >= 2) 
        Vx[0][0] += Intervald(m_a[i]*m_x[i]*(m_x[i]-1))*
                    (x[0].pow(m_x[i]-2)*x[1].pow(m_y[i])  *x[2].pow(m_z[i]));
      if (m_y[i] >= 2) 
        Vx[1][1] += Intervald(m_a[i]*m_y[i]*(m_y[i]-1))*
                    (x[0].pow(m_x[i])  *x[1].pow(m_y[i]-2)*x[2].pow(m_z[i]));
      if (m_z[i] >= 2) 
        Vx[2][2] += Intervald(m_a[i]*m_z[i]*(m_z[i]-1))*
                    (x[0].pow(m_x[i])  *x[1].pow(m_y[i])  *x[2].pow(m_z[i]-2));
      if (m_x[i] >= 1 && m_y[i] >= 1) 
        Vx[1][0] += Intervald(m_a[i]*m_x[i]*m_y[i])*
                    (x[0].pow(m_x[i]-1)*x[1].pow(m_y[i]-1)*x[2].pow(m_z[i]));
      if (m_y[i] >= 1 && m_z[i] >= 1) 
        Vx[2][1] += Intervald(m_a[i]*m_y[i]*m_z[i])*
                    (x[0].pow(m_x[i])  *x[1].pow(m_y[i]-1)*x[2].pow(m_z[i]-1));
      if (m_x[i] >= 1 && m_z[i] >= 1) 
        Vx[2][0] += Intervald(m_a[i]*m_x[i]*m_z[i])*
                    (x[0].pow(m_x[i]-1)*x[1].pow(m_y[i])  *x[2].pow(m_z[i]-1));
    }

  // Make symmetric
  Vx[0][1] = Vx[1][0];
  Vx[1][2] = Vx[2][1];
  Vx[0][2] = Vx[2][0];

  return Vx;
}

/** Returns scalar factor of the x^i y^j z^k term.
 * \param i exponent for x.
 * \param j exponent for y.
 * \param k exponent for z.
 * \return scalar factor of the x^i y^j z^k term.
 */
double Algebraic::getCoef(int i, int j, int k)
{
  return m_a[getCoefIndex(i, j, k)];
}

/** Sets coefficient of x^i y^j z^k to a.
 * \param i exponent for x.
 * \param j exponent for y.
 * \param k exponent for z.
 * \param a Coefficient for x^i y^j z^k.
 */
void Algebraic::setCoef(int i, int j, int k, double a)
{
  m_a[getCoefIndex(i, j, k)] = a;
}

int Algebraic::getNumCoef()
{
  return m_numCoef;
}

int Algebraic::degree()
{
  return m_d;
}

/** Compute the number of terms for a degree d trivariate polynomial.
 * Degree d algebraic has (d^3 + 6d^2 + 11d + 6)/6 coefficients.
 *
 * 0: 1
 * 1: 3 + 1 = 4
 * 2: 6 + 4 = 10
 * 3: 10 + 10 = 20
 *
 * # of terms of degree d = d+1 + # of terms of degree d-1.
 * Why? Multiply x times terms of degree d-1, then all that
 * is left are terms containing only y and z. There are d+1
 * of these. E.g. yyy yyz yzz zzz
 *
 * Hence terms(d) = d+1 + terms(d-1) and terms(0) = 1.
 * This yields terms(d) = sum[i = 1..d+1] i
 * = (d+2)*(d+1)/2  (Recall sum[i = 1..n] i = (n+1)*n/2)
 * = 1/2 d^2 + 3/2 d + 1.
 *
 * The # of terms of degree d or less is thus
 *   sum[i = 0..d] terms(i)
 * = sum[i = 0..d] 1/2 d^2 + 3/2 d + 1
 * = 1/2 * d*(d+1)*(2d+1)/6 + 3/2 * d*(d+1)/2 + d+1
 *            (Recall sum[i = 1..n] i^2 = n(n+1)(2n+1)/6)
 * = 1/12 * (2d^3 + 3d^2 + d) + 3/4 * (d^2 + d) + d + 1
 * = d^3/6 + d^2 + 1 10/12 d + 1.
 */
int Algebraic::coefficients(int d)
{
  return (6 + d * (11 + d * (6 + d))) / 6;
}

/** Set the number of coefficients in m_numCoef
 * based on the current degree held in m_d.
 */
void Algebraic::calcNumCoef()
{
  m_numCoef = coefficients(m_d);
}

/** Returns coefficient-order index for x^i y^j z^k term. Resulting
 * index used for m_a, m_x, m_y, m_z, m_xPow, m_yPow and m_zPow.
 *
 * Coefficient order is constant first, then x, then y, then z. For a
 * quadric, coefficient order is: 1, x, y, x^2, xy, y^2, z, xz, yz, z^2.
 *
 * \todo: Coefficient order should be by degree, so for a quadric it should
 * be: 1 x y z x^2 xy y^2 xz yz z^2. This way changing the degree means reducing
 * or extending the parameter list without necessarily reordering the lower
 * degree elements.
 *
 * 1 \\
 *
 * x y \\
 * z
 *
 * xx xy yy \\
 * xz yz \\
 * zz
 *
 * xxx xxy xyy yyy \\
 * xxz xyz yyz \\
 * xzz yzz \\
 * zzz
 *
 * Coordinates:
 *   level: d = i+j+k (degree)
 *   horizontal: j
 *   vertical: k
 *
 * Answer is:
 *   coefficients(d) - sum[kk=1..d+1-k](kk) + j
 * = coefficients(d) - sum[kk=1..i+j+1](kk) + j (since d = i+j+k)
 * = coefficients(d) - (i+j+1)(i+j+2)/2 + j
 * = coefficients(d) - (i*i + 2*i*j + j*j + 3*i + 3*j + 2)/2 + j
 *
 * \param i exponent of x
 * \param j exponent of y
 * \param k exponent of z
 * \return coefficient-order index for x^i y^j z^k.
 */
int Algebraic::getCoefIndex(int i, int j, int k)
{
    /* return (3*j + 3*(j*j + 2*i*j + i*(i + 3)) + k*(k*k - 3*(m_d + 2)*k + 3*m_d*(m_d + 4) + 11))/6; */
  // swapped j <-> i to get correct xyz order of coefficients.
    /* return (3*i + 3*(i*i + 2*j*i + j*(j + 3)) + k*(k*k - 3*(m_d + 2)*k + 3*m_d*(m_d + 4) + 11))/6; */

  /* Using above derivation */
  return coefficients(i+j+k) - (i*i + 2*i*j + j*j + 3*i + 3*j + 2)/2 + j;
}

void Algebraic::pretty(std::string &s)
{
	std::stringstream sbuf;
	int d,i,j,k,qi;
	NameVector xyz;
	Implicit::getqname(xyz);
	if (fabs(m_a[0]) != 1.0)
		xyz[0] = "";	/* instead of "1" */
	bool notfirst = false;

	for (d = degree(); d >= 0; d--) {
		for (k = 0; k <= d; k++) {
			for (j = 0; j <= d - k; j++) {
				i = d - (j + k);
				qi = getCoefIndex(i,j,k);
				if (m_a[qi] != 0.0) {
					if (notfirst) {
						if (m_a[qi] > 0.0)
							sbuf << " + ";
						else
							sbuf << " - ";
						if (fabs(m_a[qi]) != 1.0)
							sbuf << fabs(m_a[qi]);
						sbuf << xyz[qi];
					} else {
						notfirst = true;
						if (m_a[qi] == -1.0)
							sbuf << "-";
						else if (m_a[qi] != 1.0)
							sbuf << m_a[qi];
						sbuf << xyz[qi];
					}
				}
			}
		}
	}

	s = sbuf.str();
}

/** Initialize the m_x,m_y and m_z arrays to their corresponding
 * exponents in coefficient order.
 */
void Algebraic::initPowerArrays()
{
  int i, j, k;
  int coefIndex;

  for (i=0; i<=m_d; i++) 
    for (j=0; j<=m_d-i; j++) 
      for (k=0; k<=m_d-(i+j); k++) 
        {
          coefIndex = getCoefIndex(i, j, k);
          m_x[coefIndex] = i;
          m_y[coefIndex] = j;
          m_z[coefIndex] = k;
        }
}

/** Fills m_[xyz]Pow cache with powers of x, y and z.
 * Efficiently computes powers by using results of lower powers.
 * \param v The input x,y,z.
 * \note Remembers last v to avoid recomputation.
 * \todo Could memoize (hash) v to possibly avoid more recomputation.
 * Need to build a general datastructure to handle this. Perhaps this
 * should be handled at a higher level (e.g. as in Bloomenthal's marching
 * cubes implementation).
 */
void Algebraic::initXYZPowerArrays(const gmVector3 & v)
{
  int i;

  // Check to see if m_[xyz]Pow[] arrays already contain correct answer.
  if ((m_xPow[1] == v[0]) && (m_yPow[1] == v[1]) && (m_zPow[1] == v[2])) 
    return;

  m_xPow[0] = m_yPow[0] = m_zPow[0] = 1;
  for (i=1; i<=m_d; i++) 
    {
      m_xPow[i] = m_xPow[i-1] * v[0];
      m_yPow[i] = m_yPow[i-1] * v[1];
      m_zPow[i] = m_zPow[i-1] * v[2];
    }
}

double Algebraic::dx(const gmVector3 & v)
{
  double deriv = 0.0;

  initXYZPowerArrays(v);
  for (int i=0; i<m_numCoef; i++) 
    {
      if (m_x[i] >= 1) 
        deriv += m_a[i] * m_x[i] * 
                 m_xPow[m_x[i]-1] * m_yPow[m_y[i]] * m_zPow[m_z[i]];
    }
  return deriv;
}

double Algebraic::dy(const gmVector3 & v)
{
  double deriv = 0.0;

  initXYZPowerArrays(v);
  for (int i=0; i<m_numCoef; i++) 
    {
      if (m_y[i] >= 1)
        deriv += m_a[i] * m_y[i] * 
                 m_xPow[m_x[i]] * m_yPow[m_y[i]-1] * m_zPow[m_z[i]];
    }
  return deriv;
}

double Algebraic::dz(const gmVector3 & v)
{
  double deriv = 0.0;

  initXYZPowerArrays(v);
  for (int i=0; i<m_numCoef; i++) 
    {
      if (m_z[i] >= 1) 
        deriv += m_a[i] * m_z[i] * 
                 m_xPow[m_x[i]] * m_yPow[m_y[i]] * m_zPow[m_z[i]-1];
    }
  return deriv;
}

double Algebraic::dx2(const gmVector3 & v)
{
  double deriv = 0.0;

  initXYZPowerArrays(v);
  for (int i=0; i<m_numCoef; i++) 
    {
      if (m_x[i] >= 2) 
        deriv += m_a[i] * m_x[i] * (m_x[i]-1) * 
                 m_xPow[m_x[i]-2] * m_yPow[m_y[i]] * m_zPow[m_z[i]];
    }
  return deriv;
}

double Algebraic::dy2(const gmVector3 & v)
{
  double deriv = 0.0;

  initXYZPowerArrays(v);
  for (int i=0; i<m_numCoef; i++) 
    {
      if (m_y[i] >= 2) 
        deriv += m_a[i] * m_y[i] * (m_y[i]-1) * 
                 m_xPow[m_x[i]] * m_yPow[m_y[i]-2] * m_zPow[m_z[i]];
    }
  return deriv;
}

double Algebraic::dz2(const gmVector3 & v)
{
  double deriv = 0.0;

  initXYZPowerArrays(v);
  for (int i=0; i<m_numCoef; i++) 
    {
      if (m_z[i] >= 2) 
        deriv += m_a[i] * m_z[i] * (m_z[i]-1) * 
                 m_xPow[m_x[i]] * m_yPow[m_y[i]] * m_zPow[m_z[i]-2];
    }
  return deriv;
}

double Algebraic::dxdy(const gmVector3 & v)
{
  double deriv = 0.0;

  initXYZPowerArrays(v);
  for (int i=0; i<m_numCoef; i++) 
    {
      if ((m_x[i] >= 1) && (m_y[i] >= 1)) 
        deriv += m_a[i] * m_x[i] * m_y[i] * 
                 m_xPow[m_x[i]-1] * m_yPow[m_y[i]-1] * m_zPow[m_z[i]];
    }
  return deriv;
}

double Algebraic::dxdz(const gmVector3 & v)
{
  double deriv = 0.0;

  initXYZPowerArrays(v);
  for (int i=0; i<m_numCoef; i++) 
    {
      if ((m_x[i] >= 1) && (m_z[i] >= 1)) 
        deriv += m_a[i] * m_x[i] * m_z[i] * 
                 m_xPow[m_x[i]-1] * m_yPow[m_y[i]] * m_zPow[m_z[i]-1];
    }
  return deriv;
}

double Algebraic::dydz(const gmVector3 & v)
{
  double deriv = 0.0;

  initXYZPowerArrays(v);
  for (int i=0; i<m_numCoef; i++) 
    {
      if ((m_y[i] >= 1) && (m_z[i] >= 1)) 
        deriv += m_a[i] * m_y[i] * m_z[i] * 
                 m_xPow[m_x[i]] * m_yPow[m_y[i]-1] * m_zPow[m_z[i]-1];
    }
  return deriv;
}

/** Computes dproc()/dq for an algebraic.
 *  Each term of the algebraic is of the form a x^i y^j z^k so
 *  its derivative with respect to a is just x^i y^j z^k.
 */
void Algebraic::procq(const gmVector3 & x, double* dfdq) 
{
  int i;

  initXYZPowerArrays(x);
  for (i=0; i<m_numCoef; i++) 
    dfdq[i] = m_xPow[m_x[i]] * m_yPow[m_y[i]] * m_zPow[m_z[i]];  
    // dfdq is just the x^i y^j z^k of term with the given q as coefficient
}

/** Computes dproc()/dq in Interval form for an algebraic.
 *  Each term of the algebraic is of the form a x^i y^j z^k so
 *  its derivative with respect to a is just x^i y^j z^k.
 */
void Algebraic::procq(const Box<double>& x, Intervald* dfdq)
{
  int i;

  for (i=0; i<m_numCoef; i++) 
    dfdq[i] = x[0].pow(m_x[i]) * x[1].pow(m_y[i]) * x[2].pow(m_z[i]);
}

/** Computes dgrad()/dq = d^2proc/dxdq and stores it in
 * three arrays.
 *
 * \param dfdxdq returns Intervald x components of dgrad()/dq
 * \param dfdydq returns Intervald y components of dgrad()/dq
 * \param dfdzdq returns Intervald z components of dgrad()/dq
 *
 * Arrays should be preallocated to length qlen().
 *
 * Assume f(x,y,z) = a x^i y^j z^k.
 * Then df/dx = a i x^(i-1) y^j z^k.
 * Then d^2f/dxda = i x^(i-1) y^j z^k.
 */

void Algebraic::gradq(const Box<double>& x, 
                      Intervald* dfdxdq, Intervald* dfdydq, Intervald* dfdzdq)
{
  for (int i = 0; i < m_numCoef; i++) 
    {
      dfdxdq[i] = (Intervald(m_x[i]) *
                   x[0].pow(m_x[i]-1)*x[1].pow(m_y[i])  *x[2].pow(m_z[i]));
      dfdydq[i] = (Intervald(m_y[i]) *
                   x[0].pow(m_x[i])  *x[1].pow(m_y[i]-1)*x[2].pow(m_z[i]));
      dfdzdq[i] = (Intervald(m_z[i]) *
                   x[0].pow(m_x[i])  *x[1].pow(m_y[i])  *x[2].pow(m_z[i]-1));
    }
}

void Algebraic::getq(double* q) 
{
  int i;

  for (i=0; i<m_numCoef; i++)
    q[i] = m_a[i];
}

void Algebraic::_setq(double* q)
{
  for (int i=0; i<m_numCoef; i++) 
    m_a[i] = q[i];
}

unsigned int Algebraic::qlen() 
{
  return m_numCoef;
}

void Algebraic::getqname(char** qn)
{
  int i;
  char s[32];

  for (i = 0; i < m_numCoef; i++)
    {
      qn[i] = new char[32];
      qn[i][0] = '\0';
      if (m_x[i]) 
        {
          if (m_x[i] == 1) 
            strcat(qn[i],"x");
          else 
            {
              sprintf(s,"x^%d",m_x[i]);
              strcat(qn[i],s);
            }
        }
      if (m_y[i]) 
        {
          if (!qn[i][0])
            strcat(qn[i]," ");
          if (m_y[i] == 1) 
            strcat(qn[i],"y");
          else 
            {
              sprintf(s,"y^%d",m_y[i]);
              strcat(qn[i],s);
            }
        }
      if (m_z[i]) 
        {
          if (!qn[i][0])
            strcat(qn[i]," ");
          if (m_z[i] == 1) 
            strcat(qn[i],"z");
          else 
            {
              sprintf(s,"z^%d",m_z[i]);
              strcat(qn[i],s);
            }
        }
      if (!qn[i][0])
        strcat(qn[i],"1");
    }
}

