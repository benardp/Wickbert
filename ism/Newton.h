/** 
 * This file contains the class declaration for a Newton object.  This
 * class solves N equations in N variables, over a given bounding box, and
 * uses interval arithmetic to find all solutions within that bounding box.
 *
 * @author Barton Stander and John C. Hart
 * @author Terry Fleury
 * @date March 1996 (Stander), July 2001 (Hart), March 2002 (Fleury)
 */

#ifndef NEWTON_H
#define NEWTON_H

/*********************** includes ***************************/
#include "Surface/Interval.h"
#include "Surface/Box.h"
#include "tnt/vec.h"
#include "tnt/cmat.h"
#include "Surface/IMatrix.h"
#include "Surface/Implicit/Implicit.h"
#include <stack>
#include <list>

/// Enumerates the various return results from the INewton method.
typedef enum{ NO_ROOTS=1, UNIQUE_ROOT, KEEP_TRYING, NO_PROGRESS, BAD_INVERSE }
	NewtonResult;

/// Enumerates the singularity types
typedef enum{SADDLE=1,EXTREMUM,DEGENERATE,UNCLASSIFIED} SingularityType;

/// Vector of doubles - Can be of arbitrary dimension.
typedef TNT::Vector<double> TNTPoint;

/// List of points, Typedef'ed for convenience.
typedef std::vector<TNTPoint> PointList;

/// Point List with classification
typedef std::pair<TNTPoint, SingularityType> ClassifiedPoint;
typedef std::list<ClassifiedPoint> ClassifiedPointList;

/** The main Newton algorithm class.
 *
 * This module searches an N-dimensional Box for roots of an equation.  It
 * does this by first employing a simple subdivision to get rid of regions
 * of space that definitely don't contain roots.  It then subdivides those
 * regions that might contain roots into small enough pieces to be handled
 * by a Newton's Method solver.  So there are two parts of the search
 * algorithm that you must concern yourself with: subdivision and equation
 * solving.
 * 
 * This module solves N equations in N variables, over a given bounding box,
 * and is guaranteed to find all solutions within that bounding box because
 * it uses interval arithmetic.  Basically, it finds all of the roots of a
 * set of equations.
 *
 * In matrix notation, this class solves A*(Y - Xc) - b = 0 where A is an
 * an NxN IMatrix (ie. a matrix of intervals), Xc is an N-dimensional
 * vector (of doubles), and b and Y are N-dimensional Boxes (ie. vectors of
 * intervals).  Xc represents an initial guess for the root, ie. the center
 * of the interval used in Newton's method.  Y represents the variable we
 * are solving for since we are trying to find a Box containing a root.  The
 * tricky part is that YOU must provide A and b.  These values are set in
 * the NewtonEquation() method.  You will have to override this method if
 * you want to do something other than the default here (which is nothing).
 *
 * HOW TO USE THIS CLASS AND WHAT METHODS TO OVERRIDE
 *
 * (1) NewtonSubdivisionBreak() is called on a given Box region during the
 * simple subdivision phase.  This method should return true if you know for
 * sure that there is NOT a root in that region, false if there might be a
 * root in the region.  Since the default NewtonSubdivisionBreak() method
 * always returns false, the root finding won't stop until the Box size has
 * reached fWidthTolerance.
 *
 * (2) NewtonEquation() is called during the Newton's Method phase to set
 * the variables A and b in the equation A(Y-Xc) = b.  So if you want to use
 * Newton's method to refine the root finding once the Box size has reached
 * fNewtonTolerance, be sure to override this method and return true.
 *
 * (3) NewtonBreak() is called on the Box in ALL phases of the root finding
 * process to see if we can stop searching for roots in the Box.  The
 * default method performs two checks:
 *   1: Have we found the maximum number of roots allowed?
 *   2: Is the box small enough to say that two roots in the box can be
 *      considered as a single root, and if so does this box already contain
 *      a root?
 * If either of these tests is true, then we return true to indicate that we
 * want to stop searching this box for more roots.  You may want to override
 * this behavior if you have additional information on the Box in question
 * (for example, using the gradient can we determine if there is definitely
 * no root in the Box in question).
 *
 * (4) NewtonOutput() is called every time we find a root.  This method in
 * turn calls NewtonPostProcess() which allows you to take action on the
 * root in question prior to actually saving it.  You may want to override
 * NewtonOutput() if you want to store the roots someplace other than
 * the default location fRoots, a list of Points.
 * 
 * (5) NewtonPostProcess() is called by NewtonOutput() and allows you to
 * take action on the root in question prior to actually saving it.  You
 * could override this method if you have more conditions that the root must
 * meet if it will be stored.  The default method prevents duplicates roots
 * from being output.
 *
 * (6) FindZeros() is the method you call to actually find roots over an
 * initial domain. Simple subdivision is used down to a resolution of
 * fNewtonTolerance, and then Newton's method takes over, which provides
 * quadratic convergence for most cases.
 */
class Newton
{
public:
  /// How wide the divisions can be before quitting.
  static const double STANDARD_WIDTH_TOLERANCE;
  /// Maximum number of roots we can find - prevents problems with degenerates
  static const int STANDARD_MAX_ROOTS;
  // The above three consts defined at the top of Newton.cpp

  /// Maximum error of final roots (width of the bounding box).
  double fWidthTolerance; 

  /// Two roots this close together are considered the same.
  double fSameTolerance; 
  
  /// The maximum number of roots we are allowed to save.
  unsigned int fMaxRoots;

  /** Flag for FindZeros to use Newton iteration.
   * If false, then FindZeros uses simple subdivision for the entire root
   * finding process.
   */
  int fUseNewton; 

  /** FindZeros uses subdivision on large intervals, then hones in using
   * Newton when the interval is smaller than fNewtonTolerance.
   */
  double fNewtonTolerance;

  PointList fRoots;           ///< List of roots saved by NewtonOutput().
  
  PointList::iterator plIter; ///< Convenient iterator over PointList.
 
  /// Default constructor.
  Newton();               

  /// Virtual destructor.
  virtual ~Newton(){};  
  
  /// Front end for finding zeros.
  void FindZeros(Box<double> &IX);
  
protected:
  /** During the subdivsion phase, this method should return true if you
   * know for sure that there isn't a root in this box and you want to
   * exclude it from further searching during the subdivision phase.
   * @param X The Box to search for roots.
   * @return True if we should stop searching for roots in this box.
   *         False if we should search the Box for roots.
   */
  virtual bool NewtonSubdivisionBreak(Box<double> & X)
  {
    return false;
  }

  /** Setup the system of equations to solve for during the Newton's Method
   * phase.  X is the input N-dimensional Box that you want to use when
   * finding roots.  Newton's method phase.  You must set values for the
   * parameters A and b.  Also, you must return true if you want Newton
   * solver to actually do anything.  If you return false, then it assumes
   * that you haven't set up A and b and insteads uses simple subdivision
   * for the entire process.
   * @param X The input Box (N-dimensional interval vector) to search for
   *          roots.
   * @param A The interval matrix to be returned for solving the equation.
   * @param b The interval vector to be returned for solving the equation.
   * @return True if you have set A and b and want to use Newtons' method.
   *         False to indicate no processing is necessary.
   */
  virtual bool NewtonEquation(Box<double> &X, IMatrix &A, Box<double> &b)
  {
    return false;
  }

  /**
   * This method is called just before we call NewtonOutput and allows you
   * to do extra stuff to the root (in the Box) before it gets output.  If
   * you return false from here, NewtonOutput will NOT be called, in effect
   * not saving the root.  So if you want to save the root, be sure to
   * return true.  The default method checks to see if the box in question
   * has already been output.  If so, then we return false since we don't
   * want it added again.
   * @param X The Box containing the root to process just before we call
   *          NewtonOutput().
   * @return True if the root (Box) is okay and we should call
   *         NewtonOutput().  False if you want to throw the Box away.
   */
  virtual bool NewtonPostProcess(Box<double> &X)
  {
    bool retval = true;
    if (SameRoot(X))
      return false;

    return retval;
  }

  /** 
   * Outputs a root found during the search.  The input point (actually a
   * Box which gets converted to a point using the Box::center()) is checked
   * against the current list of roots (fRoots).  If it is not found in the
   * list, then it is added to the list.
   * @param X A Box which is converted to a point for adding to the fRoot
   *          list.
   * @note This method defaults to adding the root to the internal solution
   *       list. Override it if you want to maintain this list yourself.
   */
  virtual void NewtonOutput(Box<double> & X)
  {
    AddRoot(X);
	std::cout << "Number of roots = " << fRoots.size() << std::endl;
  }

  /** User definable test to determine if root search in current interval
   * can be terminated.  The default method performs the following checks:
   *   @li Have we found the maximum number of roots allowed?
   *   @li Is the box small enough to say that two roots in the box can be
   *       considered as a single root (change fSameTolerance to set this
   *       size), and if so does this box already contain a root?
   * If either of these tests is true, then we return true to indicate
   * that we want to stop searching this box for more roots. You may want to
   * add more conditions to this method (by overriding it) if you have more
   * information on the function/Box.
   * @param X Current Box being searched.
   * @return True if the Box cannot possibly contain any roots and thus we
   * should stop looking any further.  False if the Box should be searched
   * for roots.
   */
  virtual bool NewtonBreak(Box<double> &X)
  {
    // First, check to see if we have accumulated too many roots
    if (fRoots.size() > fMaxRoots)
      return true;  // Too many roots - stop looking for more

    // Next, check the size of the box for determining "sameness" of a root
    if (X.width() > fSameTolerance)
      return false;  // Box is still big - don't quit looking for roots.
    if (SameRoot(X)) // Box is pretty small - check for "same" root
      return true;   // We already found a root in this box - stop looking
    return false;    // Otherwise, Don't quit looking here
  };

  /// Determins if root already found before.
  bool SameRoot(Box<double> & IX);

  /// Add a root to the found root list.
  void AddRoot(Box<double> IX);

private:
  void FindZerosDivide(Box<double> & IX);

  void NewtonDriver();

  NewtonResult INewton(Box<double> & IX);

  int SolveSeidel(Box<double> & X, TNTPoint Xc, IMatrix & M,
    Box<double> & b, int *unique);

  int GaussSeidelRow(Box<double> & X, TNTPoint & Xc, IMatrix & M,
    Box<double> & b, int row, int *unique);
        
  
  std::stack< Box<double> > boxes; /// Stack of boxes used when root finding requires subdivision.
};

#endif 

