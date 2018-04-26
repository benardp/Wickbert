/**
 * @file IMatrix.h
 * Routines for Matrices of Intervals.
 * @author John Hart
 * @author Terry Fleury
 */

#ifndef _IMATRIX_H
#define _IMATRIX_H

#include "Box.h"
#include "tnt/cmat.h"

/** 
 * Matrix of intervals.  This class is very similar to Box except we now
 * have a 2D matrix of intervals.  See the Box class for more information.
 */
class IMatrix : public TNT::Matrix< Interval<double> > // Note space between >'s
{
  public:
    /// Constructors
    IMatrix();
    IMatrix(int m, int n);
    IMatrix(int m, int n, Interval<double> x);
    IMatrix(TNT::Matrix<double> m);
    IMatrix(TNT::Matrix< Interval<double> > m);

    /// Math operators
    IMatrix& operator *=(const Interval<double> &c);
    IMatrix& operator /=(const Interval<double> &c);
    IMatrix& operator +=(const Interval<double> &c);
    friend IMatrix operator *(const IMatrix &m, const Interval<double> &x);
    friend IMatrix operator *(const Interval<double> &x, const IMatrix &m);
    friend Box<double> operator *(const IMatrix &m1, const Box<double> &m2);
    friend IMatrix operator *(const IMatrix &m1, const IMatrix &m2);
    friend IMatrix operator /(const IMatrix &m, const Interval<double> &x);
    friend IMatrix operator +(const IMatrix &m, const Interval<double> &x);
    friend IMatrix operator +(const Interval<double> &x, const IMatrix &m);
   
    /// Size operators
    TNT::Matrix<double> center() const;

    /// Basic output
    void print();
    friend std::ostream & operator<<(std::ostream& os, const IMatrix& m);
   
    ///Identity matrix
    void MakeIdentity();
};




/**
 * This is the 3d / Double specific version of an IMatrix.  Use this most
 * often so you don't have to declare the size of the IMatrix all the time.
 */
class IMatrix3d : public IMatrix
{
  public:
    /// Constructors
    IMatrix3d();
    IMatrix3d(TNT::Matrix< Interval<double> > m);
    IMatrix3d(gmMatrix3 mindim, gmMatrix3 maxdim);
    IMatrix3d(gmMatrix3 x);
    IMatrix3d(Interval<double> xx, Interval<double> xy, Interval<double> xz,
              Interval<double> yx, Interval<double> yy, Interval<double> yz,
              Interval<double> zx, Interval<double> zy, Interval<double> zz);
    IMatrix3d(Interval<double> x);
    IMatrix3d(double d);
    IMatrix3d(IMatrix x);
   
    /// Static member functions
    static IMatrix3d identity();
};

#endif

