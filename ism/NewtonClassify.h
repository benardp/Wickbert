/**
 * This file contains the class declaration for a NewtonClassify object.
 * @file   NewtonClassify.h
 * @date   October 28, 2001
 * @author Terry Fleury (tfleury@uiuc.edu)
 */

#ifndef NEWTONCLASSIFY_H
#define NEWTONCLASSIFY_H

#include <vector>
#include "Newton.h"

/**
 * The type of critical point based on the number of negative eigen values
 * at each point.  A MAX critical point has three decreasing eigen vectors
 * ( = 3 negative eigen values ) and a MIN critical point has three
 * increasing eigen vectors ( = 3 positive eigen values ).  Note that these
 * are labelled 'backwards' from John C. Hart's presentation on critical
 * points.
 */
typedef enum { CP_UNKNOWN=-1, CP_MAX=0, CP_SADDLE2, CP_SADDLE1, CP_MIN } CPKind;

/**
 * The NewtonClassify class is a high-level class to classify critical 
 * points.  There are two main methods in this class: classify(gmVector3)
 * and classify(PointList).  The first form find the eigenvectors, 
 * eigenvalues, and "kind" of a single point.  The second form finds
 * the eigenvectors, eigenvalues and "kinds" of a list of points.  Each
 * method stores the results in internal members.
 */
class NewtonClassify
{
  private:
    /// Determine the number of positive roots and the differentEigenIndex.
    int countPosEigenValues();

  protected:
    /// Classify critical points for this Implicit surface.
    Implicit* surface;  

  public:
    /** @name Members set by Classify(PointList)
     * @{
     */
    /// A list of the kind of each critical point in the fRoots list.
    std::vector<CPKind> fKinds;
    /// A list of the eigen vectors of each critical point in fRoots.
    std::vector<gmMatrix3> fVectors;
    /// A list of the eigen values of each critical point in fRoots.
    std::vector<gmVector3> fValues;
    /** @} end Members set by Classify(PointList) */

    /** @name Members set by Classify(gmVector3)
     * @{
     */
    /// The resulting eigen vectors after a call to classify(gmVector3).
    gmMatrix3 eigenVectors;
    /// The resulting eigen values after a call to classify(gmVector3).
    gmVector3 eigenValues;
    /// The resulting kind of root after a call to classify(gmVector3).
    CPKind criticalKind;

    /** The index of the eigen value/vector that has a different sign from
     *  the others, or -1 if all eigen values are the same sign.
     */
    int differentEigenIndex;
    /** @} end Members set by Classify(gmVector3) */

    /** @name Constructors
     * @{
     */
    /// Default constructor.
    NewtonClassify();
	/// virtual destructor.
    virtual ~NewtonClassify(){	};
    /// Constructs an object containing critical points for a given Implicit.
    NewtonClassify(Implicit*);
    /** @} end of constructors */

    /** @name Setting/Getting the Surface
     * @{
     */
    /// Changes the surface to search.
    virtual void setSurface(Implicit*);
    /// Returns the current surface being searched.
    virtual Implicit* getSurface();
    /** @} end Setting/Getting the Surface */
    
    /** @name Classifying
     * @{
     */
    /// Finds the eigen vectors/values for all critical points in PointList.
    void classify(PointList);
    /// Finds the eigen vectors and eigen values of a given point.
    void classify(gmVector3);
    /** @} end classifying */
};

#endif

