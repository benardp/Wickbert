/**
* Declaration of the light position attribute.
* @file LightPosition.h
* @date 27 May. 2005
* @author Matei N. Stroila
*/

#ifndef LIGHTPOSITION_H
#define LIGHTPOSITION_H

#include "ParticleAttribute.h"

/*! \class LightPosition
 *  \brief LightPosition is an attribute which represents a particle system's light LIGHT0 position.
 */

class LightPosition : public ParticleAttribute
{
public:
	  MAKE_PARTICLESTUFF_NAME();
	  LightPosition(Particles *ps=NULL, const std::string& name=std::string("LightPosition"));
  
	  ~LightPosition();

	  virtual void prepare();

	  /*! This method returns the current light LIGHT0 position.
     *  \return mLightPosition
     */
    void const getLightPosition(gmVector4&);

    /*! This method sets the current light LIGHT0 position.
	 * and thus also the lightposition that we return
     *  \param c Sets mLightPosition
     */
    void setLightPosition(void);
	
private:

	 /*! The current light LIGHT0 position in eye coordinates.
     */
	gmVector4 _lightPositionInput;
	bool _lightIsInWorldCoords;
	gmVector4 _lightPosition;

	//temp variable to interact with GL
	float _lightPositionf[4];

}; // end class LightPosition

//set GLLight and the actual light that we return
class LightPositionSet : public PSParamButton::Callback
{
public:
	LightPosition *lightPosAttr;
	LightPositionSet(LightPosition *me) {lightPosAttr = me;}
	virtual void onbuttonpress() {lightPosAttr->setLightPosition();}
};

#endif
