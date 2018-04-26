/** 
 * Routines for handling vectors of intervals.
 * @file Box.cpp
 * @author John Hart
 * @author Terry Fleury
 */

#include "Box.h"


/****** DEFINE THE FUNCTION BODIES FOR Box3d ******/

/// Default constructor
Box3d::Box3d() : 
Box<double>(3) 
{ 
}

/** 
 * Copy constructor from the parent class.
 * This is apparently necessary since the parent class is
 * more general and for example may contain more (or less)
 * than three elements.
 */
Box3d::Box3d(TNT::Vector< Interval<double> > A) :
Box<double>(3) 
{
  (*this)[0] = A[0];
  (*this)[1] = A[1];
  (*this)[2] = A[2];
}

/**
 * Constructs a Box3d from two gmVector3s as its opposite corners.
 */
Box3d::Box3d(const gmVector3& mincorner, gmVector3 maxcorner) : 
Box<double>(3) 
{
  for (int i = 0; i < 3; i++)  
    if (mincorner[i] <= maxcorner[i])
      (*this)[i].setInterval(mincorner[i],maxcorner[i]);
    else
      (*this)[i].setInterval(maxcorner[i],mincorner[i]);
}

/**
 * This is a convenience method to make a Box3d by using a gmVector3.
 */
Box3d::Box3d(const gmVector3 & x):
Box<double>(3)
{ 
  (*this) = Box3d(x,x);
}

/**
 * This is a convenience method to make a Box3d by using three
 * Interval<double>s as the dimensions.
 */
Box3d::Box3d(Interval<double> x, Interval<double> y, Interval<double> z) : 
Box<double>(3)
{
  (*this)[0] = x;
  (*this)[1] = y;
  (*this)[2] = z;
}

/**
 * Constructs a Box3d using an Interval as all three dimensions.
 */
Box3d::Box3d(Interval<double> v) : 
Box<double>(3)
{
  (*this) = Box3d(v,v,v);
}

/**
 * Constructs a Box3d using a double as all three dimensions.
 */
Box3d::Box3d(double d) : 
Box<double>(3)
{
  (*this) = Box3d(Interval<double>(d));
}

/**
 * Constructs a Box3d from a more general Box<double>.
 * @note If the size of x is less than 3, there will be trouble!
 */
Box3d::Box3d(const Box<double>& x) :
Box<double>(3)
{
  for (int i = 0; i < 3; i++)
    (*this)[i] = x[i];
}

/**
 * Returns the 'lower corner' of the box.
 * @return The 'lower' corner of the box as a gmVector3.
 */
gmVector3 Box3d::low()
{
  return gmVector3((*this)[0].low(),
                   (*this)[1].low(),
                   (*this)[2].low());
}

/**
 * Returns the 'upper corner' of the box.
 * @return The 'upper' corner of the box as a gmVector3.
 */
gmVector3 Box3d::high()
{
  return gmVector3((*this)[0].high(),
                   (*this)[1].high(),
                   (*this)[2].high());
}




/****** DEFINE THE FUNCTION BODIES FOR Box4d ******/

/// Default constructor
Box4d::Box4d() : 
Box<double>(4) 
{ 
}

/** Copy constructor from the parent class.
 * This is apparently necessary since the parent class is
 * more general and for example may contain more (or less)
 * than four elements.
 */
Box4d::Box4d(TNT::Vector< Interval<double> > A) : 
Box<double>(4) 
{
  (*this)[0] = A[0];
  (*this)[1] = A[1];
  (*this)[2] = A[2];
  (*this)[3] = A[3];
}

/**
 * Constructs a Box4d from two gmVector4s as opposite corners.
 */
Box4d::Box4d(gmVector4 mincorner, gmVector4 maxcorner) : 
Box<double>(4) 
{
  for (int i = 0; i < 4; i++)  
    {
      if (mincorner[i] <= maxcorner[i])
        (*this)[i].setInterval(mincorner[i],maxcorner[i]);
      else
        (*this)[i].setInterval(maxcorner[i],mincorner[i]);
    }
}

/**
 * This is a convenience method to make a Box4d by using a gmVector4.
 */
Box4d::Box4d(gmVector4 x) :
Box<double>(4)
{ 
  (*this) = Box4d(x,x); 
}

/**
 * Constructs a Box4d by using four Interval<double>s as the dimensions.
 */
Box4d::Box4d(Interval<double> x, Interval<double> y, 
             Interval<double> z, Interval<double> t) : 
Box<double>(4)
{
  (*this)[0] = x;
  (*this)[1] = y;
  (*this)[2] = z;
  (*this)[3] = t;
}

/** 
 * Constructs a Box4d using an Interval as all four dimensions.
 */
Box4d::Box4d(Interval<double> v) : 
Box<double>(4)
{
  (*this) = Box4d(v,v,v,v);
}

/** 
 * Constructs a Box4d using a double as all four dimensions.
 */
Box4d::Box4d(double d) : 
Box<double>(4)
{
  (*this) = Box4d(Interval<double>(d));
}


/** 
 * Constructs a Box4d from a more general Box<double>.
 * @note If the size of x is less than 4, there will be trouble!
 */
Box4d::Box4d(const Box<double>& x) :
Box<double>(4)
{
  for (int i = 0; i < 4; i++)
    (*this)[i] = x[i];
}

/**
 * Returns the 'lower corner' of the box.
 * @return The 'lower' corner of the box as a gmVector4.
 */
gmVector4 Box4d::low()
{
  return gmVector4((*this)[0].low(),
                   (*this)[1].low(),
                   (*this)[2].low(),
                   (*this)[3].low());
}

/** 
 * Returns the 'upper corner' of the box.
 * @return The 'upper' corner of the box as a gmVector4.
 */
gmVector4 Box4d::high()
{
  return gmVector4((*this)[0].high(),
                   (*this)[1].high(),
                   (*this)[2].high(),
                   (*this)[3].high());
}

