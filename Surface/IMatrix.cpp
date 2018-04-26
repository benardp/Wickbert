/**
 * @file IMatrix.cpp
 * Routines for Matrices of Intervals
 * @author John Hart
 * @author Terry Fleury
 */

#include "IMatrix.h"

/****** DEFINE THE FUNCTION BODIES FOR IMatrix ******/

/**
 * Default constructor.  Useful only in parameter specifications.
 */
IMatrix::IMatrix() { };

/**
 * Construct an IMatrix of size m x n with initial values 0.0.
 */
IMatrix::IMatrix(int m, int n) : 
TNT::Matrix< Interval<double> >(m,n,0.0) { }

/**
 * Construct an IMatrix of size m x n with initial values x.
 */
IMatrix::IMatrix(int m, int n, Interval<double> x) : 
TNT::Matrix< Interval<double> >(m,n,x) { }

/**
 * Construct an IMatrix from a TNT::Matrix of doubles.
 */
IMatrix::IMatrix(TNT::Matrix<double> m)
{
  initialize(m.num_rows(),m.num_cols());

  for (int i = 1; i <= num_rows(); i++)
    for (int j = 1; j <= num_cols(); j++)
      (*this)(i,j) = m(i,j);
}

/**
 * Construct an IMatrix from a TNT::Matrix of Interval<double>s.
 */
IMatrix::IMatrix(TNT::Matrix< Interval<double> > m)
{
  initialize(m.num_rows(),m.num_cols());

  for (int i = 1; i <= num_rows(); i++)
    for (int j = 1; j <= num_cols(); j++)
      (*this)(i,j) = m(i,j);
}

/**
 * Return the center of each element as a TNT::Matrix of doubles.
 * @returns A TNT Matrix of doubles where each element is the 'center' of
 * each element of the IMatrix.
 */
TNT::Matrix<double> IMatrix::center() const
{
  TNT::Matrix<double> ret(num_rows(),num_cols());

  for (int i = 1; i <= num_rows(); i++)
    for (int j = 1; j <= num_cols(); j++)
      ret(i,j) = (*this)(i,j).center();

  return ret;
}

/**
 * Multiply each element of the IMatrix by an interval.
 */
IMatrix& IMatrix::operator *=(const Interval<double> &c)
{
  for (int i = 0; i < num_rows(); i++)
    for (int j = 0; j < num_cols(); j++)
      (*this)[i][j] *= c;
  return (*this);
}

/**
 * Divide each element of the IMatrix by an interval.
 */
IMatrix& IMatrix::operator /=(const Interval<double> &c)
{
  for (int i = 0; i < num_rows(); i++)
    for (int j = 0; j < num_cols(); j++)
      (*this)[i][j] /= c;
  return (*this);
}

/**
 * Add an Interval to each element of the IMatrix.
 */
IMatrix& IMatrix::operator +=(const Interval<double> &c)
{
  for (int i = 0; i < num_rows(); i++)
    for (int j = 0; j < num_cols(); j++)
      (*this)[i][j] += c;
  return (*this);
}

/**
 * Basic print out of the IMatrix.  Appends a <CR/LF> to the end.
 */
void IMatrix::print()
{
  std::cout << this << std::endl;
}

/****** DEFINE THE friend FUNCTIONS ******/

/**
 * Standard output stream to print out the coords of the box.  Note that
 * this does not output a trailing <CR/LF>.
 */
std::ostream& operator<<(std::ostream& os, const IMatrix& m)
{
  for (int i = 1; i <= m.num_rows(); i++)
    {
      os << "[ ";
      for (int j = 1; j <= m.num_cols(); j++)
        {
          os << m(i,j);
          os << ((j < m.num_cols()) ? ", " : "]");
        }
      os << " ]";
    }
  return os;
}

/// Multiply an IMatrix by an interval.
IMatrix operator *(const IMatrix &m, const Interval<double> &x)
{
  IMatrix z = m;
  return (z *= x);
}

/// Multiply an Interval by an IMatrix.
IMatrix operator *(const Interval<double> &x, const IMatrix &m)
{
  IMatrix z = m;
  return (z *= x);
}

/// Multiply an IMatrix by a Box.
Box<double> operator *(const IMatrix &m1, const Box<double> &m2)
{
  Box<double> m(m2.size(),0);
  for(int i = 0; i < m1.num_rows(); i++)
  {
	for(int j = 0; j < m1.num_cols(); j++)
  	{
		m[i] += m1[i][j] * m2[j];
	}  
  }
  return m;
}

/// Multiply an IMatrix by a IMatrix.
IMatrix operator *(const IMatrix &m1, const IMatrix &m2)
{
  IMatrix m(m1.num_rows(),m2.num_cols());
  for(int i = 0; i < m1.num_rows(); i++)
  {
	for(int j = 0; j < m2.num_cols(); j++)
  	{
	 	for(int k = 0; k < m2.num_cols(); k++)
  		{	
			m[i][j] += m1[i][k] * m2[k][j];
		}  
  	}
 }
  return m;
}

/// Divide an IMatrix by an interval.
IMatrix operator /(const IMatrix &m, const Interval<double> &x)
{
  IMatrix z = m;
  return (z /= x);
}

/// Add an IMatrix to an interval.
IMatrix operator +(const IMatrix &m, const Interval<double> &x)
{
  IMatrix z = m;
  return (z += x);
}

/// Add an Interval to an IMatrix.
IMatrix operator +(const Interval<double> &x, const IMatrix &m)
{
  IMatrix z = m;
  return (z += x);
}

void IMatrix::MakeIdentity()
{
 for(int i = 0; i < this->num_rows(); i++)
  {
	for(int j = 0; j < this->num_cols(); j++)
  	{	
	 if(i == j) (*this)[i][j] = 1;
	   else		
	    (*this)[i][j] = 0;
  	}
 }
}


/****** DEFINE THE FUNCTION BODIES FOR IMatrix3d ******/

/// Default constructor
IMatrix3d::IMatrix3d() : 
IMatrix(3,3) { }

/// Construct an IMatrix3d from a TNT Matrix of Interval<double>s
IMatrix3d::IMatrix3d(TNT::Matrix< Interval<double> > m) 
{
  initialize(3,3);
  (*this)[0][0] = m[0][0];
  (*this)[0][1] = m[0][1];
  (*this)[0][2] = m[0][2];
  (*this)[1][0] = m[1][0];
  (*this)[1][1] = m[1][1];
  (*this)[1][2] = m[1][2];
  (*this)[2][0] = m[2][0];
  (*this)[2][1] = m[2][1];
  (*this)[2][2] = m[2][2];
}

/// Construct an IMatrix3d using two gmMatrix3s as the opposite corners
IMatrix3d::IMatrix3d(gmMatrix3 mindim, gmMatrix3 maxdim)
{
  initialize(3,3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      {
        if (mindim[i][j] <= maxdim[i][j])
          (*this)[i][j] = Interval<double>(mindim[i][j],maxdim[i][j]);
        else
          (*this)[i][j] = Interval<double>(maxdim[i][j],mindim[i][j]);
      }
}

/// Contruct an IMatrix3d from a gmMatrix3
IMatrix3d::IMatrix3d(gmMatrix3 x) 
{ 
  (*this) = IMatrix3d(x,x); 
}

/// Construct the IMatrix3d the hard way - one Interval element at a time.
IMatrix3d::IMatrix3d(Interval<double>xx, Interval<double>xy, Interval<double>xz,
                     Interval<double>yx, Interval<double>yy, Interval<double>yz,
                     Interval<double>zx, Interval<double>zy, Interval<double>zz)
{
  initialize(3,3);
  (*this)[0][0] = xx;
  (*this)[0][1] = xy;
  (*this)[0][2] = xz;
  (*this)[1][0] = yx;
  (*this)[1][1] = yy;
  (*this)[1][2] = yz;
  (*this)[2][0] = zx;
  (*this)[2][1] = zy;
  (*this)[2][2] = zz;
}

/// Construct an IMatrix3d using a single Interval as value of all elements
IMatrix3d::IMatrix3d(Interval<double> x)
{
  (*this) = IMatrix3d(x,x,x,x,x,x,x,x,x);
}

/// Construct an IMatrix3d using a double as the value of all elements
IMatrix3d::IMatrix3d(double d)
{
  (*this) = IMatrix3d(Interval<double>(d));
}

/** 
 * Construct an IMatrix3d from a more general IMatrix.
 * @note If the size of x is less than 3x3, there will be trouble!
 */
IMatrix3d::IMatrix3d(IMatrix x)
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*this)[i][j] = x[i][j];
}

/// Returns the identity matrix as an IMatrix3d
IMatrix3d IMatrix3d::identity()
{
  IMatrix3d m = IMatrix3d(0.0);
  for (int i = 0; i < 3; i++)
    m[i][i] = Interval<double>(1.0,1.0);
  return m;
}

