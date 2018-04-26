/*! \file SuggestiveContourStabilization.h

 *  \brief Defines the interface to a particle behavior class which stabilizes suggestive contours.

 */



#ifndef SUGGESTIVE_CONTOUR_STABILIZATION_H

#define SUGGESTIVE_CONTOUR_STABILIZATION_H



#include "Particles.h"
#include "ParticleBehavior.h"



class ImplicitInterrogator;

class ViewDependence;



/*! \class SuggestiveContourStabilization

 *  \brief This class stabilizes suggestive contours by pushing particles away from surface locations

 *         which are likely to cause instability.  These locations are near where the view vector and

 *         surface normal become near parallel.  In other words, these are anti-silhouettes.  This

 *         behavior takes a single parameter, a threshold angle defining the boundary between stability

 *         and instability.  The force applied is gaussian weighted; in fact, the threshold is set at

 *         two standard deviations away from parallel.

 *  \see http://www.cs.rutgers.edu/~decarlo/contour.html

 *  \see SuggestiveContourAdhesion

 */

class SuggestiveContourStabilization : public ParticleBehavior

{

  public:



    MAKE_PARTICLESTUFF_NAME();

    SuggestiveContourStabilization(Particles *ps=NULL, const std::string& name=std::string("SuggestiveContourStabilization"), double threshold = 15, double scaleStep = 1);

    /// parameters
    int qlen();
    void getq(double *q);
    void setq(double *q);
    void qname(char **qn);

    void attachAttributes();

    /** Apply the suggestive contour stabilization constraint.
     */
    void applyForce();



    /*! This method returns mThreshold, the viewpoint stability threshold, in degrees.
     *  \return mThreshold
     */
    inline double getThreshold(void) const
    {
      return mThreshold;
    } // end SuggestiveContourStabilization::getThreshold()

    /*! This method sets mThreshold.
     *  \param t Sets mThreshold
     */
    inline void setThreshold(const double t)
    {
      mThreshold = t;
    } // end SuggestiveContourStabilization::setThreshold()



    /*! This method sets mScaleStep.

     *  \param s Sets mScaleStep.

     */

    inline void setScaleStep(const double s)

    {

      mScaleStep = s;

    } // end SuggestiveContourStabilization::setScaleStep()



    /*! This method returns mScaleStep.

     *  \param s Returns mScaleStep.

     */

    inline double getScaleStep(void) const

    {

      return mScaleStep;

    } // end SuggestiveContourStabilization::getScaleStep()



  protected:

    /// SuggestiveContourStabilization needs access to the surface

    ImplicitInterrogator *mImplicitInterrogator;



    /// SuggestiveContourStabilization needs access to the view

    ViewDependence *mViewDependence;
	
	ParticlePosition *position;
	ParticleVelocity *velocity;


    /// This is the threshold between stable and unstable.  Specified in degrees.

    double mThreshold;



    /// This scales the step size

    double mScaleStep;

}; // end class SuggestiveContourStabilization



#endif // SUGGESTIVE_CONTOUR_STABILIZATION_H