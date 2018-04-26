/**
 * Declaration of a contour shader 
 * @file ShaderContour.h
 * @date 09/08/2005
 * @author Matei N. Stroila
 * @remarks
 */

#ifndef SHADERCONTOUR_H
#define SHADERCONTOUR_H

#include "libgm/gmVec3.h"
#include "Particles.h"
#include "ParticleShader.h"
#include "ism/Newton.h"
#include "ism/ContourCriticalPointsGSL.h"
#include "ism/ContourCritical.h"

#include "Attributes/Contours.h"

class ParticlePosition;
class ParticleVelocity;
class ImplicitInterrogator;
class ViewDependence;
class ParticleLocality;
class AdaptiveRepulsionData;
class ParticleVector;
class ParticleBoundingBox;
class ParticleVisibility;
class Implicit;

class ShaderContour : public ParticleShader
{
public:
	enum LoopMethod {SIMPLE_LOOPS, HEAP_LOOPS};

protected:
	
	//required attributes
	ParticlePosition *position;
	ParticleVector *tangentAttr;
	ParticleVelocity *velocity;
	ViewDependence *view;
	ImplicitInterrogator *impInt;
	ParticleVisibility *vis;
	ParticleLocality *p_locality;
	AdaptiveRepulsionData *ardata;
	ParticleBoundingBox *pBoundingBox;
	Contours *contours;
		
	
	//internal variables
	bool useLoops;
	LoopMethod loopMethod;
	bool cameraPosChanged, updateConnectivity, critPointFlag, useLoopsCGAL, findCPoints;
	bool drawLines;//flag for drawing lines
	
	gmVector3 *myCameraPosition, *oldCameraPosition;
		
	unsigned int nLoops; //number of contour loops
	
	int numberCP; //number of contour critical points
	ClassifiedPointList cpList; //list of critical points
	//euler angles and zoom 
	//mouse left button down
	double xrot_last; 
	double yrot_last;
	double zoom_last;
	//mouse left button up
	double xrot; 
	double yrot;
	double zoom;
	
	
	double lineWidth; //contour width (for opengl)
	gmVector4 color, cpColor, invColor; //contours and critical points color and invisible contour colors
	//SVG variables
	gmVector3 brushColorf; //brush color
	int penWidth;
	int invPenWidth;
	bool isShaded;
	
	Particles* parabolicContoursPar; //reference to the particle collection that adheres to the parabolic contours
	std::string parabContoursName;	///< name of the parabolic contours Particles
	
	double status; //used for the crit point finding status
	
	bool useParabolicCurvesSolver, use4DSolver, use5DSolver, use_kGkR_4DSolver, useIntervalSolver; //what solver to use for finding critical points
	bool useNewton; //use the Newton method or just subdivision
	
	double NewtonTolerance, singThreshold;
	
	//internal methods

	virtual void findTangent(const unsigned int i){}; //find the tangent vector to the silhouette at the ith particle
	virtual void checkCriticalPoint4D(double *pstatus){};
	virtual void checkCriticalPoint5D(double *pstatus){};
	virtual void checkCriticalPointPC(){}; //check for critical points using the parabolic contours 
	virtual void checkCriticalPoint_kGkR_4D(double *pstatus){}; //check for crit point using the Gaussian and the radial curvature
	virtual void checkCriticalPointInterval(){}; //check for crit point using a pure interval method
	virtual void event(int);
	
	void printCPs();
	
public:
		
	ShaderContour(Particles *ps, const std::string& name);
	~ShaderContour();
	// find chains (loops)
	void findChains();
	
	virtual void drawPre();
	virtual void drawPost();

	/** Called after particles are added to the system.
	*/ 
	virtual void particleAdded();

	/** Called after particle are removed from the system.
	*/
	virtual void particleRemoved(unsigned int i);
	
	const std::vector<Loop>& getLoops() const
	{
		return contours->loops;
	}
	
	bool useVisibility;
	
	void ClearCPList(){
		cpList.clear();
	}
	
	void getContourColor(std::vector<short>& contourColor) const{
		contourColor.push_back((short) (color[0] * 255.));
		contourColor.push_back((short) (color[1] * 255.));
		contourColor.push_back((short) (color[2] * 255.));
	}

	void getInvisibleContourColor(std::vector<short>& invContourColor) const{
		invContourColor.push_back((short) (invColor[0] * 255.));
		invContourColor.push_back((short) (invColor[1] * 255.));
		invContourColor.push_back((short) (invColor[2] * 255.));
	}
	
	void getBrushColor(std::vector<short>& brushColor) const{
		brushColor.push_back((short) (brushColorf[0] * 255.));
		brushColor.push_back((short) (brushColorf[1] * 255.));
		brushColor.push_back((short) (brushColorf[2] * 255.));
	}
	const short getPenWidth() const{
		return (short) penWidth;
	}
	const short getInvPenWidth() const{
		return (short) invPenWidth;
	}
	const bool IsShaded() const{
		return isShaded;
	}
	
private:
	
	typedef std::pair<unsigned int, double> IntDoublePair;	
	struct less_dist : public std::binary_function<IntDoublePair, IntDoublePair, bool> {
		bool operator()(IntDoublePair x, IntDoublePair y) { return x.second < y.second; }
    };

	class LoopCallback:public PSParamComboBox::Callback
	{
		ShaderContour::LoopMethod * _lm;
	public: 
		LoopCallback(ShaderContour::LoopMethod * lm);
		virtual void itemselected(unsigned int selection);
	};
};

class ShaderContourClearCPList : public PSParamButton::Callback
{
public:
	ShaderContour *shc;
	ShaderContourClearCPList(ShaderContour *me) {shc = me;}
	virtual void onbuttonpress() {shc->ClearCPList();}
};




#endif

