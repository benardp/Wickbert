/**
 * This file contains the method definitions for a SearchCritical object.
 * @file   SearchCritical.cpp
 * @date   October 28, 2001
 * @author Terry Fleury (tfleury@uiuc.edu)
 */

#include "SearchCritical.h"

/**
 * Default constructor.  This method isn't very useful since it sets the
 * surface to be NULL which means that there is no surface for which to find
 * critical points.  This is just here so that you can declare a
 * SearchCritical value without having to initialize it at the same time.
 */
SearchCritical::SearchCritical()
{
  setSurface(NULL);
}

/**
 * Constructs an object containing critical points for a given Implicit.
 * You MUST pass in an Implicit surface.  This method assumes that you are
 * going to use the standard 3D Newton solver when finding critical points.
 * @param surf The implicit surface to be used when finding critical points.
 */
SearchCritical::SearchCritical(Implicit* surf)
{
  for (int i = 0; i < 3; i++)
    theBoxBounds[i].setInterval(-CP_CUBE_BOUNDS,CP_CUBE_BOUNDS);
  setSurface(surf);
}

/**
 * Changes the surface to search.  This method is called by the constructors
 * to set up the private surface variable.  It is also virtual and should be
 * overridden (along with getSurface()) to keep the surface variable
 * consistent.
 * @param surf The new Implicit surface to search for critical points.
 * @see getSurface()
 */
void SearchCritical::setSurface(Implicit* surf)
{
  surface = surf;
}

/**
 * Returns the current surface being searched.  This method is virtual and
 * should be overridden (along with setSurface(Implicit*)) to keep the
 * surface variable consistent.
 * @return The Implicit surface being searched for critical points.
 * @see setSurface(Implicit*)
 */
Implicit* SearchCritical::getSurface()
{
  return surface;
}

/**
 * Sets the internal value of the bounds for a cubic box.  You must
 * pass in a positive value for this method to do anything.  The value
 * you pass in the the length of a side of the box. So if you want
 * a box with mincorner = (-1,-1,-1) and maxcorner = (1,1,1), you 
 * should pass in bounds = 2.
 * @param bounds The length of a side of a cube to be used as the
 *               bounds for searching for critical points.
 */
void SearchCritical::setBoxBounds(double bounds)
{
  if (bounds > 0)
    for (int i = 0; i < 3; i++)
      theBoxBounds[i].setInterval(-bounds/2.0,bounds/2.0);
}

/**
 * Sets the internal value of the bounds for a cubic box.  You must
 * pass in a positive value for this method to do anything.
 * @param bounds The new bounds for searching for critical points.
 */
void SearchCritical::setBoxBounds(Box3d bounds)
{
  theBoxBounds = bounds;
}

/** 
 * Returns the current bounds for a cubic box for critical points.
 * @return The critical point bounds for the current Implicit.
 */
Box3d SearchCritical::getBoxBounds()
{
  return theBoxBounds;
}

/**
 * Overridden NewtonPostProcess method to try to prevent critical sets.
 * This is pretty much a hack, but maybe it will help in the long run.  We
 * check to see if we have reached the maximum number of roots we want to
 * store.  If so, we first increase the size of the interval used for
 * checking to see if two roots are the same, and then we double the number
 * of allowable roots.  
 * @param X The Box containing the root to process just before we call
 *          NewtonOutput().
 * @return True if the root (Box) is okay and we should call
 *         NewtonOutput().  False if you want to throw the Box away.
 */
bool SearchCritical::NewtonPostProcess(Box<double> &X)
{
  bool retval = true;
  if (SameRoot(X))
    retval = false;

  std::cout << "Num roots=" << fRoots.size()+1 << std::endl;
  // Check to see if we have already reached the maximum number of roots.
  if (fRoots.size() >= fMaxRoots)
    {
      // We probably have a critical set somewhere.  So, increase the
      // maximum number of allowable roots AND the interval size which
      // determines if two roots are the same.  But don't go overboard
      // because we don't want too many roots or have them too sparse.
      if (fRoots.size() < (fMaxRoots * 16))
        {
          if (savedfMaxRoots == 0)    // Save fMaxRoots and fSameTolerance
            {                         // for restoration after root finding
              savedfMaxRoots = fMaxRoots;
              savedfSameTolerance = fSameTolerance;
            }
          fMaxRoots *= 2;
          fSameTolerance *= 10;
          std::cout << "NEW MAX ROOTS = " << fMaxRoots << std::endl;
        }
    }
    
  return retval;
}

/**
 * Finds all of the critical points of the given surface.  Use this method
 * to do the actual work of finding the critical points of the Implicit
 * surface that you passed in at construction time.  The resulting critical
 * points will be stored in the fRoots list (found in Newton.h).  If you
 * change the parameters of your Implicit surface, you will have to
 * completely redo the finding process.
 */
void SearchCritical::search()
{
  if (getSurface() == NULL)
    return;

  savedfMaxRoots = 0;
  fNewtonTolerance = 0.1;
  Box<double> searchBox = getBoxBounds();
  FindZeros(searchBox);

  // Restore the original values of fMaxRoots and fSameTolerance
  if (savedfMaxRoots > 0)
    {
      fMaxRoots = savedfMaxRoots;
      fSameTolerance = savedfSameTolerance;
      savedfMaxRoots = 0;
      savedfSameTolerance = 0;
    }
}

/**
 * Classifies the critical points (eigen vectors/values).  This is
 * simply a convenience method to call classify(fRoots).  Make sure
 * that you have actually FOUND the roots before trying to classify
 * them.
 * @see search()
 */
void SearchCritical::classify()
{
  if (getSurface() == NULL)
    return;

  NewtonClassify::classify(fRoots);
}

/**
 * Prints all critical points to stdout.  The first line lists the name of
 * the Implicit object, and subsequent lines list the coordinates of the
 * critical points.
 */
void SearchCritical::print()
{
  int i;
  int dim;

  if (getSurface() == NULL)
    return;

  std::cout << "Critical Points for " << getSurface()->getObjectName() 
            << ":" << std::endl;
  for (plIter = fRoots.begin(); plIter != fRoots.end(); ++plIter)
    {
      std::cout << "( ";
      dim = plIter->size();
      for (i = 0; i < dim; i++)
        {
          std::cout << plIter->operator[](i);
          if (i != dim)
            std::cout << ", ";
        }
      std::cout << " )" << std::endl;
    }
}

/**
 * Saves the OpenGL state variables which are changed by display()
 * so they can be restored to their former values.
 */ 
//void SearchCritical::saveGLStateVariables()
//{
//  lightingOn = glIsEnabled(GL_LIGHTING);
//  glGetFloatv(GL_POINT_SIZE,&savedPointSize);
//  glGetFloatv(GL_CURRENT_COLOR,savedColor);
//}

/**
 * Restores the previously saved OpenGL state variables which 
 * were modified within display().
 */
//void SearchCritical::restoreGLStateVariables()
//{
//  lightingOn ? glEnable(GL_LIGHTING) : glDisable(GL_LIGHTING);
//  glPointSize(savedPointSize);
//  glColor4fv(savedColor);
//}

/**
 * Displays all critical points in an OpenGL window.
 */
//void SearchCritical::display()
//{
//  if (getSurface() == NULL)
//    return;
//
//  double sign;
//  std::vector<CPKind>::iterator kiter = fKinds.begin();
//
//  // Save all of the GL state variables
//  saveGLStateVariables();
//
//  // Now set the GL state variables for this method
//  glDisable(GL_LIGHTING);
//  glColor3f(0.9,0.9,0.9);
//
//  for (plIter = fRoots.begin(); plIter != fRoots.end(); ++plIter)
//    {
//      if (kiter!=fKinds.end())
//        switch (*kiter)
//          {
//            case CP_MAX: glColor3f(0.2,1.0,1.0); break;
//            case CP_SADDLE2: glColor3f(0.2,0.2,1.0); break;
//            case CP_SADDLE1: glColor3f(0.2,1.0,0.2); break;
//            case CP_MIN: glColor3f(1.0,0.2,0.2); break;
//          }
//      sign = getSurface()->proc(gmVector3(plIter->operator[](0),
//        plIter->operator[](1),plIter->operator[](2)));
//      glPushMatrix();
//      glTranslated(plIter->operator[](0),plIter->operator[](1),
//                   plIter->operator[](2));
//      glutSolidCube((GLdouble)((sign > 0) ? CP_POINT_SIZE : CP_POINT_SIZE/2.0));
//      glPopMatrix();
//      if (kiter!=fKinds.end())
//        ++kiter;
//    }
//
//  // Restore the previous GL state variables
//  restoreGLStateVariables();
//}

