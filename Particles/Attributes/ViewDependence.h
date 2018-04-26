/**
* Declaration of the view-dependence attribute.
* @file SilhouetteDetector.h
* @date 19 May. 2005
* @author  Jared Hoberock and Matei N. Stroila
* @brief This fine defines the interface to a ParticleAttribute encapsulating the notion of view-dependence.
*/

#ifndef VIEW_DEPENDENCE_H
#define VIEW_DEPENDENCE_H

//the camera position sensitivity
#define CAMERAPOS_SENS 0.01

#include "ParticleAttribute.h"
#include "cleangl.h"

/*! \class ViewDependence
 *  \brief ViewDependence is an attribute which represents a particle system's dependence on camera position.
 */

class ViewDependence : public ParticleAttribute
{
	friend class ViewDependenceJumpToFixCamCallback;
public:
	  MAKE_PARTICLESTUFF_NAME();

	  /// Add particle orientation to a system of particles.
	  ViewDependence(Particles *ps=NULL, const std::string& name=std::string("ViewDependence"));
	  
	  ~ViewDependence();
	  virtual void prepare();
    /*! This method returns the current camera position.
     *  \return mCameraPosition
     */

    gmVector3* getCameraPosition(void);

    /*! This method sets the current camera position.
     *  \param c Sets mCameraPosition
     */

    inline void setCameraPosition(const gmVector3 &c)
    {
      *mCameraPosition = c;
    } // end ViewDependence::setCameraPosition()

private:

class ViewDependenceJumpToFixCamCallback : public PSParamButton::Callback
{
public:
	ViewDependence *vd;
	ViewDependenceJumpToFixCamCallback (ViewDependence *me) {vd = me;}
	virtual void onbuttonpress();
};


  /*! This method computes the current camera position.
     *  \param c Computes mCameraPosition
     */

    void computeCameraPosition(void);

	 /*! The current camera position.
     */
    gmVector3* mCameraPosition;
	bool mFixCam;
	GLdouble mModelview[16];


}; // end class ViewDependence

#endif
