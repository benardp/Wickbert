/**
 * This file contains the class declaration for a SearchCritical object.
 * @file SearchCritical.h
 * @date   September 21, 2001
 * @author Terry Fleury (tfleury@uiuc.edu)
 */

#ifndef SEARCHCRITICAL_H
#define SEARCHCRITICAL_H

#define CP_CUBE_BOUNDS 1.0
#define CP_POINT_SIZE  0.05

// this is a library should not have opengl calls
// #include <GL/glut.h>
#include "NewtonCritical.h"
#include "NewtonClassify.h"

/**
 * SearchCritical is a wrapper for NewtonCritical to search for critical
 * points.  Simply call the search() method to invoke Newton's method for
 * finding roots.  The roots will be stored in the fRoots list (found in
 * Newton.h).  You can also print out the found roots to stdout by calling
 * the print() method.  Also, since SearchCritical inherits NewtonClassify,
 * you can classify all of the roots found by "search()" and get their
 * eigenvalues and eigenvectors.
 */
class SearchCritical : public virtual NewtonCritical, 
                       public virtual NewtonClassify
{
  private:
    Box3d theBoxBounds;
    //GLfloat savedPointSize;
    //GLfloat savedColor[5];
    //GLboolean lightingOn;

    /// If we change fMaxRoots, save the original value for restoration.
    int savedfMaxRoots;
    /// If we change fSameTolerance, save the original value for restoration.
    double savedfSameTolerance;

  protected:
    Implicit* surface;  ///< Find critical points for this Implicit surface.

  public:
    /// Default constructor.
    SearchCritical();
    /// Constructs an object containing critical points for a given Implicit.
    SearchCritical(Implicit*);

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
    /// Classifies the critical points (eigen vectors/values).
    void classify();
    /// Prints all critical points to stdout.
    void print();
    /// Displays all critical points in an OpenGL window.
    void display();

  protected:
    /// Overridden NewtonPostProcess method to try to prevent critical sets.
    virtual bool NewtonPostProcess(Box<double> &X);

  private:
    void saveGLStateVariables();
    void restoreGLStateVariables();
};

#endif

