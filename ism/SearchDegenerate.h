/**
 * This file contains the class declaration for a SearchDegenerate object.
 * @file SearchDegenerate.h
 * @date   September 21, 2001
 * @author Terry Fleury (tfleury@uiuc.edu)
 */

#ifndef SEARCHDEGENERATE_H
#define SEARCHDEGENERATE_H

#define CP_CUBE_BOUNDS 1.0

#include "NewtonDegenerate.h"

/**
 * SearchDegenerate is a wrapper for NewtonDegenerate to search for
 * degenerate critical points.  Simply call the search() method to invoke
 * Newton's method for finding roots.  The roots will be stored in the
 * fRoots list (found in Newton.h).  You can also print out the found roots
 * to stdout by calling the print() method.
 */
class SearchDegenerate : public NewtonDegenerate
{
  private:
    Box3d theBoxBounds;
    /// If we change fMaxRoots, save the original value for restoration.
    int savedfMaxRoots;
    /// If we change fSameTolerance, save the original value for restoration.
    double savedfSameTolerance;

  protected:
    Implicit* surface;  ///< Find critical points for this Implicit surface.

  public:
    /// Default constructor.
    SearchDegenerate();
    /// Constructs an object containing critical points for a given Implicit.
    SearchDegenerate(Implicit*);

    /// Changes the surface to search.
    virtual void setSurface(Implicit*);
    /// Returns the current surface being searched.
    virtual Implicit* getSurface();

    /// Sets the internal value of the bounds for a cubic box.
    void setBoxBounds(double);
    /// Sets the internal value of the bounds for a Box.
    void setBoxBounds(Box3d);
    /// Returns the current bounds for a Box for critical points.
    Box3d getBoxBounds();

    /// Finds all of the critical points of the given surface.
    void search();
    /// Prints all critical points to stdout.
    void print();

  protected:
    /// Overridden NewtonPostProcess method to try to prevent critical sets.
    virtual bool NewtonPostProcess(Box<double> &X);
};

#endif

