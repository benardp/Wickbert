/*
 *  ShaderSVG.h
 *  GPS
 *
 *  Created by Matei Stroila on 12/19/05.
 *  
 *
 */
#ifdef WB_USE_CGAL
#ifndef SHADERSVG_H
#define SHADERSVG_H

#include "ParticleShader.h"
#include "ParticleSystem.h"
#include "wx/wx.h"
#include "wx/dcsvg.h"
#include "boost/multi_array.hpp"
#include "PolyLines.h"
#include "ViewDependence.h"


#include <CGAL/circulator.h>

typedef  CGAL::Circulator_from_container<Loop> LoopCirculator;

typedef  CGAL::Circulator_from_container<std::vector<Point_2> > LoopPoint_2Circulator;

class ShaderSVG : public ParticleShader
{
	
public:
	
	MAKE_PARTICLESTUFF_NAME();
	
	ShaderSVG(Particles *ps=NULL);
	~ShaderSVG(){};
	
	void startCapture();
	void stopCapture();
	void capture();
	
	
	
private:
	
	typedef enum{ SILHOUETTES=1, SHADOWS, SUGGESTIVES}
	ContourType;
		
	gmVector3 penColor, brushColor, shadeBrushColor;
	int penWidth;
	
	void drawSVG(wxSVGFileDC& svgDC);
	void drawSVGContours(wxDC& dc, wxPen& pen, wxBrush& brush, std::pair<Particles* const, int>& p_pair, ContourType type, std::vector<Point_2>& pts, std::vector<Loop>& loops);
	void drawShade(wxSVGFileDC& svgDC, std::pair<Particles* const, int>& pSil_pair, std::pair<Particles* const, int>& pSha_pair);
	
	bool projectLoops(std::pair<Particles* const, int>& p_pair, ContourType type, std::vector<Point_2>& pts, std::vector<Loop>& loops);
	
	void computeVisualEvents(std::pair<Particles* const, int>& pSil_pair, std::pair<Particles* const, int>& pSha_pair);
	
	PolyLines *polyLines;
	std::string filename;
	bool captureShaded;
	
	Pmwx myPM; //planar map with intersection
	std::vector<Point_2> silPts;
	std::vector<Point_2> shadPts;
	std::vector<Point_2> suggPts;
		
	std::vector<Loop> silLoops;
	std::vector<Loop> shadLoops;
	std::vector<Loop> suggLoops;
	
	void computePlanarMap(wxSVGFileDC& dc, std::pair<Particles* const, int>& pSil_pair, std::pair<Particles* const, int>& pSha_pair);
	
	ViewDependence *viewDep;
	
	bool findClosestPoint(const Point_2& sourcePt,  Point_2& targetPt); //return true if found
	bool findClosestParticle(const Point_2& sourcePt, 
							const std::vector<Loop>&, 
							std::pair<Particles* const, int>& pSil_pair,
							Point_2& targetPt,
							int& p_index, //the closest particle index
							int& l_index); //the corresponding loop )

};

class ShaderSVGCapture : public PSParamButton::Callback
{
public:
	ShaderSVG *sha;
	ShaderSVGCapture(ShaderSVG *me) {sha = me;}
	virtual void onbuttonpress() {sha->capture();}
};


#endif
#endif