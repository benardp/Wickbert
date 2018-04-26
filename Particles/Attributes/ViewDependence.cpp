/*! \file ViewDependence.cpp
 *  \author Jared Hoberock
 *  \brief Implementation of ViewDependence class.
 */

#include "ViewDependence.h"

REGISTER_PARTICLESTUFF(ViewDependence,"Attribute:ViewDependence");

void ViewDependence::ViewDependenceJumpToFixCamCallback::onbuttonpress()
{
	glMatrixMode(GL_MODELVIEW_MATRIX);
	glLoadMatrixd(vd->mModelview);
}

ViewDependence::ViewDependence(Particles *ps, const std::string& name):ParticleAttribute(ps,name)
{
	mCameraPosition = new gmVector3(0.0,0.0,0.0);
	new PSParamgmVector3(this,mCameraPosition,gmVector3(0,0,0),
		"camera","camera position","Position of camera used for visibility computation.");
	
	new PSParamBool(this, &mFixCam, false,"fixCam","fix view","fixes the current viewpoint");
	new PSParamButton(this,new ViewDependenceJumpToFixCamCallback(this),"jumpBTN","jump to fix position","Place the camera at the fixed/current position");
} // end ViewDependence::ViewDependence()


void ViewDependence::prepare()
{
	if (!mFixCam)
		computeCameraPosition();
}

void ViewDependence::computeCameraPosition(void)
{
	glGetDoublev(GL_MODELVIEW_MATRIX,mModelview);
	gmMatrix4 _m(mModelview[0], mModelview[4], mModelview[8], mModelview[12],
		mModelview[1], mModelview[5], mModelview[9], mModelview[13],
		mModelview[2], mModelview[6], mModelview[10], mModelview[14],
		mModelview[3], mModelview[7], mModelview[11], mModelview[15]);
	gmMatrix4 _m_inv = _m.inverse();
	(*mCameraPosition)[0] = _m_inv[0][3];
	(*mCameraPosition)[1] = _m_inv[1][3];
	(*mCameraPosition)[2] = _m_inv[2][3];
}// end ViewDependence::computeCameraPosition()

gmVector3* ViewDependence::getCameraPosition(void)
{
	return mCameraPosition;
} // end ViewDependence::getCameraPosition()

ViewDependence::~ViewDependence()
{
	 delete mCameraPosition;
}