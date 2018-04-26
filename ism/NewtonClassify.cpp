/**
 * This file contains the method definitions for a NewtonClassify object.
 * @file   NewtonClassify.cpp
 * @date   October 30, 2001
 * @author Terry Fleury (tfleury@uiuc.edu)
 */

#include <math.h>
#include "NewtonClassify.h"
#include "Jacobi.h"

/**
 * Default constructor.  This method isn't very useful since it sets the
 * surface to be NULL which means that there is no surface for which to
 * classify critical points.  This is just here so that you can declare a
 * NewtonClassify value without having to initialize it at the same time.
 */
NewtonClassify::NewtonClassify()
{
  setSurface(NULL);
}

/**
 * Constructs an object containing classification of critical points for a
 * given Implicit. You MUST pass in an Implicit surface. 
 * @param surf The implicit surface to be used when classifying critical
 *             points.
 */
NewtonClassify::NewtonClassify(Implicit* surf)
{
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
void NewtonClassify::setSurface(Implicit* surf)
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
Implicit* NewtonClassify::getSurface()
{
  return surface;
}

/**
 * Finds the eigen vectors/values for all points in the roots list.  The
 * eigen vectors of the points are stored in the fVectors list.
 * The eigen values of the points are stored in the fValues list.
 * The kinds of the points are stored in the fKinds list.  
 * @param roots A list of roots previously found by Newton which we want to
 *              classify.
 * @note This method sets the object members fKinds, fVectors, and fValues.
 *       These three lists correspond to the kinds, eigenvectors, and
 *       eigenvalues of the passed-in roots.
 */
void NewtonClassify::classify(PointList roots)
{
  PointList::iterator plIter; 

  if (getSurface() == NULL)
    return;

  fKinds.clear();
  fVectors.clear();
  fValues.clear();
  for (plIter = roots.begin(); plIter != roots.end(); ++plIter)
    {
      classify(gmVector3(plIter->operator[](0),
                         plIter->operator[](1),
                         plIter->operator[](2)));
      fVectors.push_back(eigenVectors);
      fValues.push_back(eigenValues);
      fKinds.push_back(criticalKind);
    }
}

/** 
 * Finds the eigen vector and eigen value of a given point.  The resulting
 * values are stored in the class members eigenVectors and eigenValues.
 * Also, this method figures out the "kind" of the point by counting the
 * number of positive eigen values and setting the criticalKind member.
 * It also stores the index of the eigen value that had a different sign
 * from the rest, or -1 if all eigen values were of the same sign.
 * @param p The point used to find eigen vectors/values with the current
 *          Implicit surface.
 * @note This method sets the object members eigenvectors, eigenvalues,
 *       criticalKind, and differentEigenIndex.
 */
void NewtonClassify::classify(gmVector3 p)
{
  if (getSurface() == NULL)
    return;

  gmMatrix3 hess = getSurface()->hess(p);
  Jacobi jacob(hess);
  jacob.solve();
  jacob.sort();
  eigenVectors = jacob.vectors;
  eigenValues = jacob.values;
  criticalKind = (CPKind)(countPosEigenValues());
}

/** 
 * Determine the number of positive roots and the differentEigenIndex. This
 * is called by classify(gmVector3) after the eigen values are found.
 * @return The number of positive eigen values in the eigenValues member.
 * @note This method also sets the object's differentEigenIndex value which
 *       indicates which eigen value was a "different sign" from the others,
 *       or -1 if all of the eigen values are all the same sign.
 */
int NewtonClassify::countPosEigenValues()
{
  if (getSurface() == NULL)
    return -1;

  int poscount = 0;
  int negcount = 0;
  int posindex, negindex;

  for (int d = 0; d < 3; d++)
    {
      if (eigenValues[d] > 0.0)
        {
          poscount++;
          posindex = d;
        }
      else if (eigenValues[d] < 0.0)
        {
          negcount++;
          negindex = d;
        }
    }

  if (poscount == 1)
    differentEigenIndex = posindex;
  else if (negcount == 1)
    differentEigenIndex = negindex;
  else
    differentEigenIndex = -1;

  return poscount;
}

