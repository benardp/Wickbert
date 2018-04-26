/** Routines for handling vectors of intervals.
 * @file Box.h 
 * @author John Hart
 * @author Terry Fleury
 * @date 3 July 2001
 */

#ifndef BOX_H
#define BOX_H

#include "tnt/vec.h"
#include "Interval.h"
#include "RandomStream.h"

/** 
 * @class Box 
 * An interval vector class.
 * A box is an axis aligned region of space represented by a vector of
 * intervals.  A generalized Box can be of any dimension.  Specialized
 * classes exist for 3D and 4D regions.
 * 
 * @todo Any need for this to depend on TNT:vector instead of std::vector?
 */
template <class Float>
class Box : public TNT::Vector< Interval<Float> > // Note space between >'s
{
  private:
    RandomStream rs;   ///< A stream of random numbers

  public:
    /// Constructors
    Box();
    Box(int n);
    Box(int n, Float x);
    Box(TNT::Vector<Float> x);
    Box(TNT::Vector< Interval<Float> > x);

    /// Math operators
    Box<Float>& operator *=(const Interval<Float> &c);
    Box<Float>& operator /=(const Interval<Float> &c);
    Box<Float>& operator +=(const Interval<Float> &c);

    /// Size operators
    Float width(int &dim) const;  ///< Also returns which side was longest
    Float width() const;
    Float diagonal() const;
    Interval<Float> length() const;
    Interval<Float> lengthSquared() const;
    Box<Float> normalize();

    /// Parts of Boxes
    Box<Float> unionWith(Box<Float> &b);
    void subdivide(Box<Float> &IX1, Box<Float> &IX2);
    Box<Float> intersection(Box<Float> &b) const;
    void intersection(Box<Float> &b);
    TNT::Vector<Float> center() const;

    /// Basic output
    void print();
};  // End of class definition




/// Explicit instantiation of Box<double>
//template class Box<double>;

/// Multiply a box and an interval.
template <class Float>
Box<Float> operator *(const Box<Float> &x, const Interval<Float> &y)
{
  Box<Float> z = x;
  return (z *= y);
}

/// Multiply an interval and a box.
template <class Float>
Box<Float> operator *(const Interval<Float> &x, const Box<Float> &y)
{
  Box<Float> z = y;
  return (z *= x);
}

/// Divide a box and an interval.
template <class Float>
Box<Float> operator /(const Box<Float> &x, const Interval<Float> &y)
{
  Box<Float> z = x;
  return (z /= y);
}

/// Add a box and an interval.
template <class Float>
Box<Float> operator +(const Box<Float> &x, const Interval<Float> &y)
{
  Box<Float> z = x;
  return (z += y);
}

/// Add an interval and a box.
template <class Float>
Box<Float> operator +(const Interval<Float> &x, const Box<Float> &y)
{
  Box<Float> z = y;
  return (z += x);
}

/// 'Dot' product of two boxes - loose definition of dot product here.
template <class Float>
Interval<Float> dot(const Box<Float> &x, const Box<Float> &y)
{
  Interval<Float> sum = 0.0;
  for (int i = 0; i < x.size(); i++)
    sum += x[i] * y[i];

  return sum;
}

/// Standard output stream to print out the coords of the box.
template<class Float>
std::ostream& operator<<(std::ostream& os, const Box<Float>& x)
{
  os << "[ ";
  for (int i = 0; i < x.size(); i++)
    {
      os << x[i];
      if (i < (x.size()-1))
        os << ", ";
    }
  return os << " ]";
}




/**
 * This is the 3D / double specific version of a Box.  Use this most often
 * so you don't have to declare the size of your Boxes all the time.
 */
class Box3d : public Box<double>
{
  public:
    /// Constructors
    Box3d();
    Box3d(TNT::Vector< Interval<double> > A);
    Box3d(const gmVector3& mincorner, gmVector3 maxcorner);
    Box3d(const gmVector3 & x);
    Box3d(Interval<double> x, Interval<double> y, Interval<double> z);
    Box3d(Interval<double> v);
    Box3d(double d);
    Box3d(const Box<double>& x);

    /// 'Corner' operators
    gmVector3 low();
    gmVector3 high();
};




/**
 * This is the 4D / double specific version of a Box.  Use this when you
 * need a 4D box (ie, x,y,z,t) for finding 4D critical points (through time).
 */
class Box4d : public Box<double>
{
  public:
    /// Constructors
    Box4d();
    Box4d(TNT::Vector< Interval<double> >);
    Box4d(gmVector4 mincorner, gmVector4 maxcorner);
    Box4d(gmVector4 x);
    Box4d(Interval<double>,Interval<double>,Interval<double>,Interval<double>);
    Box4d(Interval<double>);
    Box4d(double d);
    Box4d(const Box<double>& );
   
    /// 'Corner' operators
    gmVector4 low();
    gmVector4 high();
};

/****** DEFINE THE FUNCTION BODIES FOR Box<Float> ******/

/**
* Default constructor.  Used for parameter specifications.
 */
template <class Float>
Box<Float>::Box() : 
TNT::Vector< Interval<Float> >() { }

/**
* Construct a Box of size n.
 * @param n The dimensionality of the Box.
 */
template <class Float>
Box<Float>::Box(int n) : 
TNT::Vector< Interval<Float> >(n) { }

/**
* Construct a Box of size n with initial values x.
 * @param n The dimensionality of the Box.
 * @param x The initial values for all dimensions of the Box.
 */
template <class Float>
Box<Float>::Box(int n, Float x) : 
TNT::Vector< Interval<Float> >(n,Interval<Float>(x)) { }

/**
* Construct a Box from a TNT vector of doubles.  This is a convenience
 * operator to convert a TNT vector to a Box.
 * @param x The TNT::Vector to transform into Box.
 */
template <class Float>
Box<Float>::Box(TNT::Vector<Float> x) : 
TNT::Vector< Interval<Float> >(x.size()) 
{
	for (int i = 0; i < x.size(); i++) 
		(*this)[i] = Interval<Float>(x[i]);
}

/**
* Construct a Box from a TNT vector of intervals.  
 * @param x The TNT::Vector< Interval<Float> > to transform into Box.
 */
template <class Float>
Box<Float>::Box(TNT::Vector< Interval<Float> >x) : 
TNT::Vector< Interval<Float> >(x.size()) 
{
	for (int i = 0; i < x.size(); i++) 
		(*this)[i] = x[i];
}

/**
* Multiply the Box times an Interval.  Each element of the box is
 * multiplied by the scalar.  Amazingly, TNT doesn't have this operator.
 * @param c An Interval to multiply with each element of this box. 
 */
template <class Float>
Box<Float>& Box<Float>::operator *=(const Interval<Float> &c)
{
	for (int i = 0; i < this->size(); i++)
		(*this)[i] *= c;
	return (*this);
}

/**
* Divide the Box by an Interval.  Each element of the box is divided by the
 * scalar.
 * @param c An Interval to divide into each element of this box. 
 */
template <class Float>
Box<Float>& Box<Float>::operator /=(const Interval<Float> &c)
{
	for (int i = 0; i < this->size(); i++)
		(*this)[i] /= c;
	return (*this);
}

/**
* Add an Interval to each element of the Box.  
 * @param c An Interval to add to each element of this box. 
 */
template <class Float>
Box<Float>& Box<Float>::operator +=(const Interval<Float> &c)
{
	for (int i = 0; i < this->size(); i++)
		(*this)[i] += c;
	return (*this);
}

/** 
* Find the width of longest side of box.  This method takes a 'dummy'
* parameter which gets set to the side (index) with the largest dimension
* upon return.
* @param dim Returns the widest dimension.
* @return The size of widest dimension.
*/
template <class Float>
Float Box<Float>::width(int &dim) const
{
	Float maxwidth = -1.0;
	
	for (int i = 0; i < this->size(); i++) 
		if ((*this)[i].width() > maxwidth) 
		{
			dim = i;
			maxwidth = (*this)[i].width();
		}
			return maxwidth;
}

/**
* Convenience function to find the maximum length of the sides of the box.
 * @see width(int &dim)
 * @return The size of the widest dimension.
 */
template <class Float>
Float Box<Float>::width() const
{
	int dummy;
	return width(dummy);
}

template <class Float>
Float Box<Float>::diagonal() const
{
    Float retVal, temp;

	for (int i = 0; i < this->size(); i++) 
    {
        temp = (*this)[i].width();
        retVal = temp*temp;
    }

    return sqrt(retVal);
}

/**
* This method is poorly named, but we need some equivalent to "length"
 * of a vector.  So given a Box such as [ [xlo,xhi],[ylo,yhi],[zlo,zhi] ]
 * we want to return an Interval which gives the closest/and furthest
 * distance from the origin to the interval.  
 * @return An interval of the [closest,furthest] points of the Box.
 */
template <class Float>
Interval<Float> Box<Float>::length() const
{
	Interval<Float> lensq = (*this).lengthSquared();
	return lensq.sqrt();
}

/**
* This method is poorly named, but we need some equivalent to
 * "lengthSquared" of a vector.  So given a Box such as 
 * [ [xlo,xhi],[ylo,yhi],[zlo,zhi] ] we want to return an Interval which
 * gives the closest/and furthest distance from the origin to the
 * interval, squared.  
 * @return An interval of the [closest,furthest]^2 points of the Box.
 */
template <class Float>
Interval<Float> Box<Float>::lengthSquared() const
{
	Infinite<Float> lo = 0;
	Infinite<Float> hi = 0;
	
	for (int i = 0; i < this->size(); i++)
    {
		if ((*this)[i].isNegative() || (*this)[i].isPositive())
			lo += (Infinite<Float>::min(fabs((*this)[i].low()),
										fabs((*this)[i].high()))).pow(2);
		hi += (Infinite<Float>::max(fabs((*this)[i].low()),
									fabs((*this)[i].high()))).pow(2);
    }
	return Interval<Float>(lo,hi);
}

/**
* We need some equivalent to "normalize" for a standard vector.  So, we
 * divide each element of the Box by it's "length" and set the box equal to
 * these new elements.
 * @return The 'normal' of the box = Box / 'length of Box'
 */
template <class Float>
Box<Float> Box<Float>::normalize()
{
	Interval<Float> len = (*this).length();
	*this /= len;
	return *this;
}

/** 
* Returns the union of two boxes.
* @param b A box to union with the current box.
* @return A new box enclosing the two original boxes.
*/
template <class Float>
Box<Float> Box<Float>::unionWith(Box<Float> &b) 
{
	Box<Float> retbox;
	for (int i = 0; i < this->size(); i++)
		retbox[i] = (*this)[i].unionWith(b[i]);
	return retbox;
}

template <class Float>
Box<Float> Box<Float>::intersection(Box<Float> &b) const
{
	Box<Float> retbox;
	for (int i = 0; i < this->size(); i++)
		retbox[i] = (*this)[i].intersection(b[i]);
	return retbox;
}

template <class Float>
void Box<Float>::intersection(Box<Float> &b)
{
	for (int i = 0; i < this->size(); i++)
		(*this)[i] = (*this)[i].intersection(b[i]);
}

/** 
* Subdivide a box nearly in half along widest dimension.  The two new
* 'halves' of the Box are returned in IX1 and IX2.
* 
* @param IX1 Left half of interval
* @param IX2 Right half of interval
*
* @note Actual subdivision happens randomly between 45% and 55% of the
* width.  A random breaking point (hopefully) avoids problems where the
* split point happens to be exactly the solution point.
*/
template <class Float>
void Box<Float>::subdivide(Box &IX1, Box &IX2)
{
	int i;
	Float w = width(i);  // Note this sets i!!!
#if 1
	Float ratio = rs.next(0.4,0.6);
#else
	Float ratio = .5;
#endif
	Float middle = (1.0 - ratio) * (*this)[i].low() + ratio * (*this)[i].high();
	
	IX1 = IX2 = *this;
	
	IX1[i] = Interval<Float>((*this)[i].low(),middle);
	IX2[i] = Interval<Float>(middle,(*this)[i].high());
	
}

/**
* Return the 'center' of the box as a TNT Vector.
 * @return A TNT::Vector of the 'center' of the Box.
 */
template <class Float>
TNT::Vector<Float> Box<Float>::center() const
{
	int i = 0;
	typename Box<Float>::iterator elem;
	TNT::Vector<Float> center(this->size());
	
	for (elem = this->begin(); elem != this->end(); elem++) 
		center[i++] = elem->center();
	
	return center;
}

/**
* Basic print out of the box.  Appends a <CR/LF> to the end.
 */
template <class Float>
void Box<Float>::print()
{
	std::cout << (*this) << std::endl;
}





#endif

