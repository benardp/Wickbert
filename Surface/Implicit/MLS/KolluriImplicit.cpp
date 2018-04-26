/**
 * @file KolluriImplicit.cpp
 *
 * @brief Implementation of KolluriImplicit class.
 */
#ifdef WB_USE_SFL

#include "KolluriImplicit.h"

using namespace MLS;

// Registers the implicit the with factory.
REGISTER_IMPLICIT(KolluriImplicit,"MLS:KolluriImplicit");

/**
 * @brief Creates an KolluriImplicit.
 *
 * Creates an KolluriImplicit.  Sets the normal and energy functions we will use for
 * the MLS.
 */
KolluriImplicit::KolluriImplicit(void) :  MLSSurface() {

    _normalF = MLSSurface::computeOptMLSDir;
    _energyF = MLSSurface::computeMLSEnergy;
 }



/**
 * @brief Computes the value of the implicit at x.
 * 
 * Uses Adamson and Alexa's definition of an Implicit form of
 * an MLS Surface to define the value of the implicit at x.
 */
double KolluriImplicit::proc(const gmVector3 & x) {

    double avgDist; //implicit value

    //compute the weighted average distance from the sample points in a
    //neighborhood around x
    if (MLSSurface::computeDistanceAvg(x, m_Kernel, avgDist))
        return avgDist;
    else
        return 0.;
}

/**
 * Returns the gradient (normal) at x.  We use the
 * averages of the Surflet normals.
 */
gmVector3 KolluriImplicit::grad(const gmVector3 & x) {
    
    Vector norm(0, 0, 0);
    Point avg;  //return values

    MLSSurface::computeAvgNormal(x, m_Kernel, norm, avg);
    
    return norm;
}

/**
 * Evaluates the gradient everywhere in the box x and returns this as a three dimension Interval (box).
 */

Box3d KolluriImplicit::grad(const Box<double>& x) {

    Box3d retVal;
    std::vector<DistancedIndex> neighbors;
  
    //The neighboring points as defined by the kernel
    m_Kernel->getNeighbors(x, m_Tree, neighbors);
        
    retVal[0] = pde1Box(x, neighbors, 0);
    retVal[1] = pde1Box(x, neighbors, 1);
    retVal[2] = pde1Box(x, neighbors, 2);

    return retVal;
}

/**
 * Evaluates the partial derivative of proc with respect to the argth coordinate over the domain x
 * given sample point neigbors.
 */

Interval<double> KolluriImplicit::normalPde1Box(const Box<double>& x,
    std::vector<DistancedIndex> neighbors, int arg) {

    Interval<double> retVal(0, 0);
    int found = neighbors.size();
  
    //Now collect all weighted normals
    Interval<double> w, totWeight(0., 0.);
    Sample s;
    for (int i = 0; i < found; ++i) {
        s = m_Samples[neighbors[i].index];
        w = (*m_Kernel)(x, s);
        totWeight += w;
        retVal += w*s.n[arg];
    }

    //to little weight, use the closest one, unweighted
    if (totWeight.high() < 0.0001) {
        std::nth_element(neighbors.begin(), neighbors.begin(), neighbors.end());
        s = m_Samples[neighbors.begin()->index];
        retVal = Interval<double>(s.n[arg], s.n[arg]);
    }

    return retVal;
}

/**
 * Evaluates the partial derivative of proc with respect to the argth coordinate over the domain x
 * given sample point neighbors.  The PDE includes the weight.
 */

Interval<double> KolluriImplicit::weightedPde1Box(const Box<double>& x,
    std::vector<DistancedIndex> neighbors, int arg) {

    Interval<double> retVal(0, 0);
    int found = neighbors.size();
    GaussianClip *kernel = dynamic_cast<GaussianClip*>(m_Kernel);

    Interval<double> totWeight(0, 0);
    Interval<double> totWeightD(0, 0);
    for (int i = 0; i < found; ++i) {
        totWeight += (*m_Kernel)(x, m_Samples[neighbors[i].index]);
        totWeightD += m_Kernel->pde1Box(x, m_Samples[neighbors[i].index], arg);
    }

    Interval<double> w;
    Sample s;
    for (int i = 0; i < found; ++i) {
        s = m_Samples[neighbors[i].index];
        w = (*m_Kernel)(x, s);

        retVal += 
            (w*s.n[arg]) / totWeight
            - 
            ((2.0*w) * (x[arg] - s[arg]) * (s.n[0]*(x[0] - s[0]) + s.n[1]*(x[1] - s[1]) + s.n[2]*(x[2] - s[2]))) /
            (kernel->getSqrVariance()*totWeight)
            -
            (w * (s.n[0]*(x[0] - s[0]) + s.n[1]*(x[1] - s[1]) + s.n[2]*(x[2] - s[2])) * totWeightD) /
            totWeight.squared();
    }

    retVal /= kernel->getSqrVariance() * totWeight;
    return retVal;
}

/**
 * Evaluates the partial derivative of proc with respect to the argth coordinate over the domain x
 * given sample point neighbors.
 */

Interval<double> KolluriImplicit::pde1Box(const Box<double>& x,
    std::vector<DistancedIndex> neighbors, int arg) {

    Interval<double> retVal(0, 0);
    int found = neighbors.size();
    GaussianClip *kernel = dynamic_cast<GaussianClip*>(m_Kernel);

    Interval<double> w, totWeight(0., 0.);
    Interval<double> t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
    Sample s;
    for (int i = 0; i < found; ++i) {

        s = m_Samples[neighbors[i].index];
        w = (*kernel)(x, s);

        t1 = w;
        t2 = s.n[arg];
        t3 = 2.0 * (x[arg] - s[arg]);
        t4 = s.n[0]*(x[0] - s[0]);
        t5 = s.n[1]*(x[1] - s[1]);
        t6 = s.n[2]*(x[2] - s[2]);
        t7 = t4 + t5 + t6;
        t8 = t3 * t7;
        t9 = t8/(kernel->getSqrVariance());
        t10 = t2 - t9;
        t11 = t1 * t10;

        totWeight += w;
        retVal += w * 
            (s.n[arg] - 
                (2.0 * (x[arg] - s[arg]) * (s.n[0]*(x[0] - s[0]) + s.n[1]*(x[1] - s[1]) + s.n[2]*(x[2] - s[2]))) /
                (kernel->getSqrVariance()));
    }
#if 0
    //to little weight, use the closest points with no weighting done to them
    if (totWeight.high() < 0.00001) {

        //find the radius of the next query, box diagonal size + radius to nearest point
        m_Tree->nearestQuery(convert(x.center()), 1, neighbors);
        double r = sqrt(neighbors.begin()->dist2) + x.diagonal()/2.;

        //find the closest set of samples with the KDTree
        m_Tree->radiusQuery(convert(x.center()), r, neighbors);

        //iterate over every neighbor and find the biggest n[arg] and smallest n[arg]
        double minN = 10, maxN = -10;   //assumed to be unitized, so they should be between 1 and -1
        for (std::vector<DistancedIndex>::iterator i = neighbors.begin(); i != neighbors.end(); ++i) {

            s = m_Samples[i->index];
            minN = s.n[arg] < minN ? s.n[arg] : minN;
            maxN = s.n[arg] > maxN ? s.n[arg] : maxN;
        }
        
        retVal = Interval<double>(minN, maxN);
    }
#endif
    return retVal;
}

#endif

