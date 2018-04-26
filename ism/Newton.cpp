/** @file Newton.cpp Interval Newton Solver
 * @author Barton Stander and John C. Hart
 * @date 30 Dec. 1995 (originally)
 * @brief Object-oriented implementation of an interval Newton system solver.
 *
 * Use FindZeros to find ALL zeros of N equations in N variables
 *
 * This module solves N equations in N unknowns, over a given bounding box,
 * and is guaranteed to find all solutions withing that bounding box because
 * it uses interval arithmetic.  Simple subdivision is used down to a
 * resolution of fNewtonTolerance, and then newton's method takes over,
 * which provides quadratic convergence in many cases.
 *
 * For documentation on how to use the methods in this class, see the
 * Newton.h header file.
 */

#include "Newton.h"

// Definition of Newton's class constants
const double Newton::STANDARD_WIDTH_TOLERANCE = 0.001;
const int Newton::STANDARD_MAX_ROOTS =20;

/** Default constructor.  It initializes the algorithm variables.
 */
Newton::Newton()
{
  fWidthTolerance = STANDARD_WIDTH_TOLERANCE;
  fSameTolerance = fWidthTolerance * 100.0;
  fMaxRoots = STANDARD_MAX_ROOTS;
  fUseNewton = true;
  // For better performance, USER should set this one!
  fNewtonTolerance = 0.5; 
}

/** Determines if root already found before.  The "sameness" is determined
 * by the member fSameTolerance which is related to fWidthTolerance.
 * @param IX Small Box interval surrounding newly found root.
 * @return True if nearby root already found, false otherwise.
 */
bool Newton::SameRoot(Box<double> & IX)
{
  TNTPoint diff;
  TNTPoint newroot = IX.center();
  bool same;

  for (plIter = fRoots.begin(); plIter != fRoots.end(); ++plIter) 
    {
      same = true;   // Assume new point is already in the list of roots
      for (int d = 0; d < IX.size(); d++)
        if (fabs(plIter->operator[](d) - newroot[d]) > fSameTolerance)
          {
            same = false;  // The new point and the current root are different 
            break;
          }
      if (same)
        return true;
    }
  return false;
}

/** Add a root to the found root list.
 * \param IX small interval around newly found root
 */
void Newton::AddRoot(Box<double> IX) 
{ 
  fRoots.push_back(IX.center()); 
}


/** Perform a step in the Gauss-Seidel iteration for a given row.
 * @param X Domain of root finding (ie. the Box used)
 * @param Xc Initial guess (ie. the midpoint of the Box)
 * @param M Interval matrix = A in the matrix equation Ax = b
 * @param b Interval vector = b in the matrix equation Ax = b
 * @param row Current row to work on.
 * @param unique Set to true if interval contains at most one root.
 *        Undefined if no solutions found.
 * @return
 *   @li False if there are no solutions here.
 *   @li True if X contains at least one box to be further processed.
 *
 * @note All values here are Interval values
 */
int Newton::GaussSeidelRow(Box<double> & X, TNTPoint & Xc, 
                           IMatrix & M, Box<double> & b, 
                           int row, int *unique)
{
  int dimen = X.size();    // Dimensionality of the Box X
  Interval<double> R;      // Intermed. result = b-SUM(Mij(Xj-xj)) j=1..N,j!=i
  Interval<double> Y1,Y2;  // Possible solutions for Ri / Mii
  Interval<double> T1,T2;  // Intersections of possible solutions
  int j,split,ok1,ok2;     // Loop iterators and logic results
  
  /* Initialize R to b_row to calculate the quantity                    */
  /* R_row = b_row - SUM (M_ij * (X_j - Xc_j)) for j=0..dimen-1, j!=row */
  R = b[row];
  for (j = 0; j < dimen; j++)
    if (j != row)
      R -= ( (X[j]-Xc[j]) * M[row][j] );

  /* Perform R_row / M_row_row using Interval division */
  /* split = IDinf< Interval<double> >(Y1,Y2,R,M[row][row]); */
  split = IDinf(Y1,Y2,R,M[row][row]); 

  /* Check first interval for possible intersection for solution */
  Y1 += Xc[row];
  if (ok1 = X[row].overlaps(Y1))
    T1 = X[row].intersection(Y1);

  /* If the Interval division resulted in a split, check second interval */
  ok2 = false;
  Y2 += Xc[row];
  if (split == ISPLIT)
    if (ok2 = X[row].overlaps(Y2))
      T2 = X[row].intersection(Y2);

  /* Handle all possible cases for ok1 and ok2 */
  if(!ok1 && !ok2)
    return false ; /* No solutions in X */

  if (!ok2) 
    { /* T2 is empty, but T1 is ready to go */
      *unique = X[row].contains(Y1);
      X[row] = T1;
      return true;
    }

  if (!ok1) 
    { /* T1 is empty, but T2 is ready to go */
      *unique = X[row].contains(Y2);
      X[row] = T2;
      return true;
    }

  /* ok1 && ok2 => Both halfs of Y intersect X */
  /* Deal with T2 sub-box later */
  X[row] = T2;
  boxes.push(X);

  /* Continue with element T1 in box X */
  X[row] = T1; 

  *unique = false;  // Both intervals half-infinite, can't be a subset of X.

  return true;
} /* End GaussSeidelRow */

/** Solve the interval matrix problem M(Y-Xc) + b = 0 for Y.
 * Set X = Intersection(X,Y).
 * Split X and push a sub-box if necessary.
 * @param X Domain of the iteration (ie. the Box used)
 * @param Xc Seed point of the iteration (ie. the center of Box)
 * @param M Interval matrix = A in the matrix equation Ax = b
 * @param b Interval vector = b in the matrix equation Ax = b
 * @param *unique True if X contains a single root
 * @return
 *   @li False if there are no solutions here.
 *   @li True if X contains a box to be further processed.
 *   @li if return==true && *unique==true, then X contains a unique root
 * @note: All values here are Interval values
 */
int Newton::SolveSeidel(Box<double> &X, TNTPoint Xc, IMatrix &M,
                        Box<double> &b, int *unique)
{
  int dimen = X.size();         // Number of dimensions in the Box
  
  int *gaps = new int[dimen];   // Find when diag matrix entries contain 0.0
  int i;                        // Loop iterator
  int sub;                      // Is the subinterval unique?

  /* Assume Y is a subinterval of X in all dimensions */  
  *unique=true;

  /* Check for when diagonal of matrix M crosses 0.0 */
  for (i = 0; i < dimen; i++)
    gaps[i] = M[i][i].contains(0.0);
  
  /* Process X in dimensions where M[i][i] does NOT cross 0 */
  for (i = 0; i < dimen; i++) 
    if (!gaps[i]) 
      {
        if (!GaussSeidelRow(X,Xc,M,b,i,&sub))
          return false; /* No solutions in X */
        if (!sub)
          *unique = false;
      }

  /* Process X in dimensions where M[i][i] DOES cross 0.0 */
  for (i = 0; i < dimen; i++)
    if (gaps[i]) 
      {
        if (!GaussSeidelRow(X,Xc,M,b,i,&sub))
          return false; /* No solutions in X */
        if (!sub)
          *unique = false;
      }

  delete [] gaps;
  /* Got this far => X needs further processing */
  return true;
} /* End SolveSeidel */

/** Interval newton method for solving simultaneous equations of the 
 * form A*(Y-XC) - b = 0.
 * @param X Box in which to find roots.
 * @return
 *   @li NO_ROOTS if there is no solution in IX
 *   @li UNIQUE_ROOT if IX contains a unique solution
 *   @li KEEP_TRYING if not sure, but IX is getting smaller
 *   @li NO_PROGRESS if not sure, but IX is not getting smaller
 *   @li BAD_INVERSE if Jc cannot be inverted
 * @see NewtonEquation() for setting up A and b in the matrix equation Ax=b.
 */
NewtonResult Newton::INewton(Box<double> &X)
{
  IMatrix A;                   // Jcinv*J (approximately identity)
  Box<double> b;               // Jcinv*Fc;
  Box<double> Xold;            // Original Box X
  int unique;                  // Uniqueness of solution

  /* User must supply the matrix A and vector b in the equation */
  /* A(Y-Xc) - b = 0 where we are solving for Y. */
  if (!NewtonEquation(X,A,b))
    return(NO_ROOTS);

  Xold = X;

  /* For A(Y-Xc) - b = 0             */
  /* we must solve for the new box Y */

  if(!SolveSeidel(X,(TNTPoint)(X.center()),A,b,&unique)) /* No solutions in X */
    return NO_ROOTS;

  /* Check for necessary subdivision */
  if(X.width() > 0.9*Xold.width()) /* not significantly smaller */
    return NO_PROGRESS;

  /* Check for a unique solution */
  if (unique) /* Y is a subinterval of X */
    return UNIQUE_ROOT;

  /* KEEP_TRYING- don't know yet, keep iterating this function */
  return KEEP_TRYING;
} /* End INewton */

/** Performs the Newton zero search on whatever is in the stack.
 */
void Newton::NewtonDriver()
{
  Box<double> X;
  int iterationsnewton;

  while (!boxes.empty()) /* more boxes to examine */
    {
      X = boxes.top();
      boxes.pop();
      iterationsnewton = 0;

      while (true) 
        {
          iterationsnewton++;

          if (iterationsnewton > 20) 
            {
              std::cout<<"Warning! More than 20 Newton Iterations"<<std::endl;
              break;
            }

          NewtonResult status=INewton(X);

#if 0  
          std::cout << "Status = ";
          switch (status) 
            {
              case NO_ROOTS: std::cout << "NO_ROOTS"; break;
              case UNIQUE_ROOT: std::cout << "UNIQUEROOT"; break;
              case KEEP_TRYING: std::cout << "KEEP_TRYING"; break;
              case NO_PROGRESS: std::cout << "NO_PROGRESS"; break;
              case BAD_INVERSE: std::cout << "BAD_INVERSE"; break;
              default: std::cout << "UNKNOWN"; break;
            }
          std::cout << std::endl;
#endif
            
          if (NewtonBreak(X)) // User wants to skip this box
            break;

          if (X.width() < fWidthTolerance) 
            { // Found close enough, quit
#if 0
              if (status != UNIQUE_ROOT) // We'll probably want to break here
                std::cout<<"Warning - found root is not unique!"<<std::endl;
#endif

              if (NewtonPostProcess(X))
                NewtonOutput(X);
              break;
            }

          if(status==NO_ROOTS) // Found nothing, quit
            break;

          if(status==NO_PROGRESS || status==BAD_INVERSE) 
            { // Subdivide the box and work on each half
              Box<double> X1;
              Box<double> X2;
              X.subdivide(X1,X2);
              boxes.push(X1);
              boxes.push(X2);
              break;
            }
        }
    }
} /* End NewtonDriver */

/** 
 * This method iterates through the given Box region looking for roots.  It
 * first uses a simple subdivision method to search for roots.  If
 * fUseNewton is set, it switches over to Newton's Method when the Box is
 * sufficiently small (set by fNewtonTolerance).
 * @param X Box region to search for roots (zeros).
 */
void Newton::FindZerosDivide(Box<double> &X)
{
  Box<double> X1,X2; /* The two halves if we divide X */
  double width;      /* The width of the input box */

  /* Check if root is definitely NOT found in the current Box */
  if(NewtonSubdivisionBreak(X))
    return;
            
  /* See if the user wants to quit this box */
  if(NewtonBreak(X))
    return;

  /* Use Newton's method on this box if it is small enough. */
  width = X.width();
  if ((fUseNewton) && (width < fNewtonTolerance)) 
    {
      boxes.push(X);
      NewtonDriver();
      return;
    }

  /* As a last resort, if the interval is small enough, call it a root! */
  if (width < fWidthTolerance) 
    {
      if (NewtonPostProcess(X))
        NewtonOutput(X);
      return;
    } 

  /* Unsure, so subdivide */
  X.subdivide(X1,X2);
  FindZerosDivide(X1);
  FindZerosDivide(X2);
}

/** Front end for finding zeros.  This is the method to call when you want
 * to look for roots in the given region Box.
 * @param X Box in which to find zeros.
 */
void Newton::FindZeros(Box<double> &X)
{
  fRoots.clear();       // Make sure the list of roots is empty.
  FindZerosDivide(X);   // Call the real worker method.
}

