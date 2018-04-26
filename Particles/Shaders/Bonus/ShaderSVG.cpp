/*
 *  ShaderSVG.cpp
 *  GPS
 *
 *  Created by Matei Stroila on 12/19/05.
 *  
 *
 */
#ifdef WB_USE_CGAL
#include "ShaderSVG.h"
#include "ShaderContour.h"
#include "ParticleVisibility.h"
#include "ParticlePosition.h"
#include "ParticleNormal.h"
#include "cleangl.h"


REGISTER_PARTICLESTUFF(ShaderSVG,"Shader:ShaderSVG");

ShaderSVG::ShaderSVG(Particles *ps)
:ParticleShader(ps,std::string("ShaderSVG"))
{
	new PSParamgmVector3(this,&penColor,gmVector3(0,0,0),
						 "penColor","Pen Color","Contours Pen Color");
	new PSParamInt(this,&penWidth,4,
				   "penWidth","Pen Width","Pen Width");
	new PSParamgmVector3(this,&brushColor,gmVector3(250,193,146),
						 "brushColor","Brush Color","Brush Color");
	new PSParamgmVector3(this,&shadeBrushColor,gmVector3(253,230,211),
						 "shadeBrushColor","Shade Brush Color","Shade Brush Color");
	new PSParamString(this,&filename,"wxModeler.svg","Filename","SVG File Name");
	new PSParamButton(this,new ShaderSVGCapture(this),"capture","Capture","Capture single frame SVG");
	new PSParamBool(this,&captureShaded,true,"Shaded","Shaded","Set false if want only contours");
	new Attached<PolyLines>(this,&polyLines);
	new Attached<ViewDependence>(this,&viewDep);
}

void ShaderSVG::capture(){

	///save the SVG to a file
	int w, h;
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	w = viewport[2];
	h = viewport[3];
	
	wxSVGFileDC svgDC (filename, w, h);
	drawSVG(svgDC);
	
}

void ShaderSVG::drawSVG(wxSVGFileDC& dc)
{
	std::map<Particles*, int> silhParticlesMap, shadowParticlesMap, suggestiveParticlesMap; //maps particle collection/index of the contour shader in the shader vector
	Particles* currentPCollection;
	int foundSilh, foundShadow, foundSugg; //index of the contour shader
	
	ParticleSystems *pss = ps->particleSystem->particleSystems;
	//find contour particle collections
	for (unsigned int j = 0; j < pss->size(); j++){
		for(unsigned int i = 0; i < (*pss)[j]->particles.size(); i++){
			currentPCollection = (*pss)[j]->particles[i];
			foundSilh = currentPCollection->findShader("ShaderSilhouetteContour");
			foundShadow = currentPCollection->findShader("ShaderShadowContour");
			foundSugg = currentPCollection->findShader("ShaderSuggestiveContour");
			if(foundSilh != -1) silhParticlesMap[currentPCollection] = foundSilh;
			if(foundShadow != -1) shadowParticlesMap[currentPCollection] = foundShadow;
			if(foundSugg != -1) suggestiveParticlesMap[currentPCollection] = foundSugg;
		}
	}
	//draw silhouettes

	std::map<Particles*, int>::iterator iter, iterSil;
	for(iter = silhParticlesMap.begin(); iter != silhParticlesMap.end(); ++iter)
	{
		wxBrush wB;
		if(captureShaded) 
			wB.SetColour(brushColor[0],brushColor[1],brushColor[2]);
		else
			wB = *wxTRANSPARENT_BRUSH;
		
		wxPen wP;
		wP.SetColour(penColor[0],penColor[1],penColor[2]);
		wP.SetWidth(penWidth);
		
		dc.SetPen(wP);
		if(!projectLoops( *(silhParticlesMap.begin()),SILHOUETTES, silPts, silLoops)) continue;
		//drawSVGContours(dc, wP, wB, *iter, SILHOUETTES, silPts, silLoops);
	}
	//draw shadows
	for(iter = shadowParticlesMap.begin(); iter != shadowParticlesMap.end(); ++iter)
	{
		wxBrush wB(*wxTRANSPARENT_BRUSH);
		if(captureShaded) 
			wB.SetColour(shadeBrushColor[0],shadeBrushColor[1],shadeBrushColor[2]);
		else
			wB = *wxTRANSPARENT_BRUSH;
		
		wxPen wP;
		wP.SetColour(penColor[0],penColor[1],penColor[2]);
		wP.SetWidth(1);
		dc.SetPen(wP);
		if(!projectLoops( *(shadowParticlesMap.begin()),SHADOWS, shadPts, shadLoops)) return;
		//drawSVGContours(dc, wP, wB, *iter, SHADOWS, shadPts, shadLoops);
	}
	
	if(captureShaded){
		drawShade(dc, *(silhParticlesMap.begin()), *(shadowParticlesMap.begin()));
	}
	//draw suggestive contours
	
	for(iter = suggestiveParticlesMap.begin(); iter != suggestiveParticlesMap.end(); ++iter)
	{
		wxPen wP;
		wP.SetColour(penColor[0],penColor[1],penColor[2]);
		wP.SetWidth(penWidth);
		dc.SetPen(wP);
		drawSVGContours(dc, wP, *wxTRANSPARENT_BRUSH, *iter, SUGGESTIVES, suggPts, suggLoops);
	}
	
}

bool ShaderSVG::projectLoops(std::pair<Particles* const, int>& p_pair, ContourType type, std::vector<Point_2>& pts,  std::vector<Loop>& loops)
{
	pts.clear();
	loops.clear();
	//get contour loops
	Particles* p_collection = p_pair.first;
	int found = p_pair.second;
	ShaderContour *sc = dynamic_cast<ShaderContour *>(p_collection->shaders[found]);
	loops = sc->getLoops();
	
	ParticlePosition *position = dynamic_cast<ParticlePosition *>(p_collection->getAttributeGeneric("ParticlePosition"));
	if(!position) return false;
	
	//if no loops present, return
	if(loops.size() == 0) return false;
	
	//project loops
	GLdouble projection[16];  
	glGetDoublev(GL_PROJECTION_MATRIX, projection);// Retrieve The Projection Matrix
	GLdouble modelview[16];  
		
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);// Retrieve The Modelview Matrix	
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	
	gmVector4 pi; //particle position
	
	//project the Silhouette loops
	int no_particles = p_collection->size();
	// Create a 2D array that is no_particles x 3
	typedef boost::multi_array<double, 2> array_type;
	array_type screenCoord(boost::extents[no_particles][3]);
	for(int i=0;i<no_particles;i++)
	{
		pi = gmVector4(position->getPosition(i)[0], position->getPosition(i)[1], position->getPosition(i)[2], 1.0);
		gluProject(
				   pi[0],
				   pi[1],
				   pi[2],
				   modelview,
				   projection,
				   viewport,
				   &(screenCoord[i][0]),
				   &(screenCoord[i][1]),
				   &(screenCoord[i][2])
				   );
		Point_2 pt(screenCoord[i][0],screenCoord[i][1], i); //corresponding screen space position
		pts.push_back(pt);
	}
	return true;
}

void ShaderSVG::drawSVGContours(wxDC& dc, wxPen& pen, wxBrush& brush, std::pair<Particles* const, int>& p_pair, ContourType type, std::vector<Point_2>& pts,  std::vector<Loop>& loops)
{
	std::vector<Loop>::const_iterator loops_iter; //iterator over loops
	Loop::const_iterator loop_iter, loop_iter_next;	//iterator over a loop specified by loops_iter
	
	//draw shaded
	dc.SetBrush(brush);
	dc.SetPen(*wxTRANSPARENT_PEN);
	int* count = new int[loops.size()]; //array of loops sizes
	int* p_count = count;
	int n_points = 0; //total number of points
	for(loops_iter = loops.begin(); loops_iter != loops.end(); ++loops_iter)
	{
		*p_count = loops_iter->size();
		n_points += *p_count;
		p_count++;
	}
	
	wxPoint* points = new wxPoint[n_points];
	
	int i = 0; //index of the current point
	for(loops_iter = loops.begin(); loops_iter != loops.end(); ++loops_iter)
	{
		for(loop_iter = loops_iter->begin(); loop_iter != loops_iter->end(); ++loop_iter)
		{
			int x0 = (int)CGAL::to_double(pts[*loop_iter].x()); 
			int y0 = (int)CGAL::to_double(pts[*loop_iter].y());
			points[i] = wxPoint(x0, y0);
			i++;
		}
	}
			
	dc.DrawPolyPolygon( loops.size(), count,  points, 0, 0, wxWINDING_RULE);
			
	//draw loops
	//get contour loops
	Particles* p_collection = p_pair.first;
	int found = p_pair.second;
	ShaderContour *sc = dynamic_cast<ShaderContour *>(p_collection->shaders[found]);
	int isVisble1 =1;
	int isVisble2 =1;
	ParticleVisibility *vis = dynamic_cast<ParticleVisibility *>(p_collection->getAttributeGeneric("ParticleVisibility"));
	if(!vis) return;
	wxPen wP = *wxTRANSPARENT_PEN; //pen for invisible parts
	dc.SetBrush(*wxTRANSPARENT_BRUSH);
	dc.SetPen(pen);
	for(loops_iter = loops.begin(); loops_iter != loops.end(); ++loops_iter)
	{
		loop_iter_next = loops_iter->begin();
		++loop_iter_next;
		for(loop_iter = loops_iter->begin(); loop_iter != loops_iter->end(); ++loop_iter)
		{
			if(loop_iter_next == loops_iter->end()) loop_iter_next = loops_iter->begin();
			
			if(sc->useVisibility){
				isVisble1 = vis->getVisibility(*loop_iter);
				isVisble2 = vis->getVisibility(*loop_iter_next);
				if((isVisble1 == 1) && (isVisble2 == 1))
				{
					dc.SetPen(pen);
				}
				else
				{
					dc.SetPen(wP);
				}
			}
			dc.DrawLine((int)CGAL::to_double(pts[*loop_iter].x()), (int)CGAL::to_double(pts[*loop_iter].y()), (int)CGAL::to_double(pts[*loop_iter_next].x()), (int)CGAL::to_double(pts[*loop_iter_next].y()));
			++loop_iter_next;
		}
	}
	delete [] points;
	delete[] count;
			
}

void ShaderSVG::computePlanarMap(wxSVGFileDC& dc, std::pair<Particles* const, int>& pSil_pair, std::pair<Particles* const, int>& pSha_pair)
{

//	Pmwx::Face_iterator fit;
//	Pmwx::Ccb_halfedge_circulator cc;
//	Pmwx::Halfedge_handle hh;
	
//	myPM.clear();
	X_monotone_curve_2::const_iterator pit, pit_next;

	dc.SetPen(*wxRED_PEN);
	std::vector<Loop>::iterator loops_iter; //iterator over loops
	Loop::iterator loop_iter, loop_iter_next; //iterator over a loop specified by loops_iter
	typedef std::vector<Point_2> Loop_pts;		//points of the loops_iter loop
	Loop_pts loop_pts;
	std::cout << "no of sil loops: " << silLoops.size() << std::endl;
	std::vector<Loop_pts> visSilhouettes, visShadows;

	int p_idx, l_idx;

	//get silhouette visibility
	Particles* p_collection = pSil_pair.first;
	int found = pSil_pair.second;
	ShaderContour *sc = dynamic_cast<ShaderContour *>(p_collection->shaders[found]);
	int isVisble1 =1;
	int isVisble2 =1;
	ParticleVisibility *vis = dynamic_cast<ParticleVisibility *>(p_collection->getAttributeGeneric("ParticleVisibility"));
	if(!vis) return;
	bool onEnd = false;
	for(loops_iter = silLoops.begin(); loops_iter != silLoops.end(); loops_iter++)
	{
		//initialize circulators
		LoopCirculator loopCc(&(*loops_iter));
		LoopCirculator loopCcNext(&(*loops_iter));
		LoopCirculator loopCcInit(&(*loops_iter));
		++loopCcNext;
		
		//find first self-intersection
		while((vis->getVisibility(*loopCcNext) == vis->getVisibility(*loopCcInit)) 
			&& (loopCcNext != loopCcInit )){
			++loopCcNext;
			++loopCc;
		}
		if(loopCcNext == loopCcInit && vis->getVisibility(*loopCcInit) == 0) continue; //no visible parts, move to the next loop 
		//if on invisible part move to the next visible part 
		if(vis->getVisibility(*loopCcNext) == 0){
			loopCcInit = loopCcNext;
			++loopCcNext;
			++loopCc;			
			while((vis->getVisibility(*loopCcNext) == vis->getVisibility(*loopCcInit)) 
				&& (loopCcNext != loopCcInit )){
				++loopCcNext;
				++loopCc;			
			}
		}
		loopCcInit = loopCc;
		do{
			onEnd = false;
			isVisble1 = vis->getVisibility(*loopCc);
			isVisble2 = vis->getVisibility(*loopCcNext);
			if((isVisble1 == 0) && (isVisble2 == 1)){
				//find closest silhouette particle
				Point_2 closestPoint;
				if(findClosestParticle(silPts[*loopCcNext], silLoops, pSil_pair, closestPoint, p_idx, l_idx))
				{
					loop_pts.push_back(	closestPoint );
					dc.SetPen(*wxBLACK_PEN);
					dc.DrawPoint((int)CGAL::to_double(closestPoint.x()), (int)CGAL::to_double(closestPoint.y()));
				}
				++loopCc;
				++loopCcNext;
			}
			if((isVisble1 == 1) && (isVisble2 == 0)){
				loop_pts.push_back(silPts[*loopCc]);
				//find closest silhouette particle
				Point_2 closestPoint;
				if(findClosestParticle(shadPts[*loopCc], silLoops, pSil_pair, closestPoint, p_idx, l_idx))
				{
					loop_pts.push_back(	closestPoint );
					dc.SetPen(*wxBLACK_PEN);
					dc.DrawPoint((int)CGAL::to_double(closestPoint.x()), (int)CGAL::to_double(closestPoint.y()));
				}
				++loopCc;
				++loopCcNext;
				onEnd = true;
			}
			while((isVisble1 == 0) && (isVisble2 == 0)){
				loopCc++;
				loopCcNext++;
				isVisble1 = vis->getVisibility(*loopCc);
				isVisble2 = vis->getVisibility(*loopCcNext);
				if(loopCc == loopCcInit) break;
			}
			while((isVisble1 == 1) && (isVisble2 == 1)){
				loop_pts.push_back(silPts[*loopCc]);
				loopCc++;
				loopCcNext++;
				isVisble1 = vis->getVisibility(*loopCc);
				isVisble2 = vis->getVisibility(*loopCcNext);
				if(loopCc == loopCcInit) {
					loop_pts.push_back(silPts[*loopCc]);
					onEnd = true;
					break;
				}
			}
			if(loop_pts.size()>1 && onEnd)
			{	
				LoopPoint_2Circulator ccLoop(&loop_pts);
				LoopPoint_2Circulator ccLoopNext(&loop_pts);
				LoopPoint_2Circulator ccInit(&loop_pts);
				ccLoopNext++;
				while(ccLoopNext != ccInit){
					dc.SetPen(*wxRED_PEN);
					dc.DrawLine((int)CGAL::to_double(ccLoop->x()), (int)CGAL::to_double(ccLoop->y()), 
						(int)CGAL::to_double(ccLoopNext->x()), (int)CGAL::to_double(ccLoopNext->y()));
					ccLoop++;
					ccLoopNext++;
				};
				visSilhouettes.push_back(loop_pts);
//				myPM.insert (Curve_2(Org_curve_2(loop_pts.begin(), loop_pts.end()),-1)); //data -1 means curves coming from silhouettes
				
			}
			if(onEnd) loop_pts.clear();
		}while(loopCc != loopCcInit);
	}
	
	//add the shadow contours loops
	std::cout << "no of shadow loops: " << shadLoops.size() << std::endl;
	//get shadow visibility
	p_collection = pSha_pair.first;
	found = pSha_pair.second;
	if(found == -1) return;
	sc = dynamic_cast<ShaderContour *>(p_collection->shaders[found]);
	isVisble1 =1;
	vis = dynamic_cast<ParticleVisibility *>(p_collection->getAttributeGeneric("ParticleVisibility"));
	if(!vis) return;
	loop_pts.clear();
	
	std::vector<Curve_2> shadowCurves; //curves to be inserted

	typedef std::pair< int, int> ShadowEnd; //loop index - particle index
	typedef std::pair<ShadowEnd, ShadowEnd> ShadowEnds; 
	ShadowEnds s_ends;
	ShadowEnd s1_end, s2_end;
	std::vector<ShadowEnds> shadowsEnds;
	std::vector<ShadowEnds>::iterator ends_iter;
	for(loops_iter = shadLoops.begin(); loops_iter != shadLoops.end(); loops_iter++)
	{
		//initialize circulators
		LoopCirculator loopCc(&(*loops_iter));
		LoopCirculator loopCcNext(&(*loops_iter));
		LoopCirculator loopCcInit(&(*loops_iter));
		++loopCcNext;			
		while((vis->getVisibility(*loopCcInit) == vis->getVisibility(*loopCcNext)) 
			&& (loopCcNext != loopCcInit )){
			++loopCcNext;
			++loopCc;			
		}
		if(loopCcNext == loopCcInit && vis->getVisibility(*loopCcInit) == 0) continue; //no visible parts, move to the next loop 
		//if on invisible part move to the next visible part 
		if(vis->getVisibility(*loopCcNext) == 0){
			loopCcInit = loopCcNext;
			++loopCcNext;
			++loopCc;			
			while((vis->getVisibility(*loopCcInit) == vis->getVisibility(*loopCcNext)) && (loopCcNext != loopCcInit )){
				++loopCcNext;
				++loopCc;			
			}
		}
		loopCcInit = loopCc;
		do{
			onEnd = false;
			isVisble1 = vis->getVisibility(*loopCc);
			isVisble2 = vis->getVisibility(*loopCcNext);
			if((isVisble1 == 0) && (isVisble2 == 1)){
				//find closest silhouette particle
				Point_2 closestPoint;
				if(findClosestParticle(shadPts[*loopCcNext], silLoops, pSil_pair, closestPoint,  p_idx, l_idx))
				{
					loop_pts.push_back(	closestPoint );
					s1_end = std::make_pair(l_idx, p_idx);
					dc.SetPen(*wxBLACK_PEN);
					dc.DrawPoint((int)CGAL::to_double(closestPoint.x()), (int)CGAL::to_double(closestPoint.y()));
				}
				++loopCc;
				++loopCcNext;
			}
			if((isVisble1 == 1) && (isVisble2 == 0)){
				loop_pts.push_back(shadPts[*loopCc]);
				//find closest silhouette particle
				Point_2 closestPoint;
				
				if(findClosestParticle(shadPts[*loopCc], silLoops, pSil_pair, closestPoint, p_idx, l_idx))
				{
					loop_pts.push_back(	closestPoint );
					s2_end = std::make_pair(l_idx, p_idx);
					dc.SetPen(*wxBLACK_PEN);
					dc.DrawPoint((int)CGAL::to_double(closestPoint.x()), (int)CGAL::to_double(closestPoint.y()));
				}
				++loopCc;
				++loopCcNext;
				onEnd = true;
			}
			while((isVisble1 == 0) && (isVisble2 == 0)){
				++loopCc;
				++loopCcNext;
				isVisble1 = vis->getVisibility(*loopCc);
				isVisble2 = vis->getVisibility(*loopCcNext);
				if(loopCc == loopCcInit) break;
			}
			while((isVisble1 == 1) && (isVisble2 == 1)){
				loop_pts.push_back(shadPts[*loopCc]);
				++loopCc;
				++loopCcNext;
				isVisble1 = vis->getVisibility(*loopCc);
				isVisble2 = vis->getVisibility(*loopCcNext);
				if(loopCc == loopCcInit) {
					loop_pts.push_back(shadPts[*loopCc]);
					s1_end = std::make_pair(-1, -1);
					s2_end = std::make_pair(-1, -1);
					onEnd = true;
					break;
				}
			}
			if(loop_pts.size()>1 && onEnd)
			{	
				shadowCurves.push_back(Curve_2(Org_curve_2(loop_pts.begin(), loop_pts.end()),0)); //data 0 means curves coming from shadows
				visShadows.push_back(loop_pts);
				s_ends = std::make_pair(s1_end, s2_end); 
				shadowsEnds.push_back(s_ends);
				LoopPoint_2Circulator ccLoop(&loop_pts);
				LoopPoint_2Circulator ccLoopNext(&loop_pts);
				LoopPoint_2Circulator ccInit(&loop_pts);
				ccLoopNext++;
				while(ccLoopNext != ccInit){
					dc.SetPen(*wxGREEN_PEN);
					dc.DrawLine((int)CGAL::to_double(ccLoop->x()), (int)CGAL::to_double(ccLoop->y()), 
						(int)CGAL::to_double(ccLoopNext->x()), (int)CGAL::to_double(ccLoopNext->y()));
					++ccLoop;
					++ccLoopNext;
				};
			}	
			if(onEnd) loop_pts.clear();
		}while(loopCc != loopCcInit);
	}

	std::vector<Loop_pts>::iterator visShadows_iter = visShadows.begin();
	Loop_pts::iterator visShadow_iter;
	for(ends_iter = shadowsEnds.begin(); ends_iter != shadowsEnds.end(); 
		ends_iter++, visShadows_iter++)
	{
		if(ends_iter->first.first == ends_iter->first.first == -1)
			continue;

		if(ends_iter->first.first == ends_iter->second.first) //same silhouette loop
		{
			int l = ends_iter->first.first;
			int p1 = ends_iter->second.second;
			int p2 = ends_iter->first.second; 

/*
			Loop::iterator iter_loop1 = silLoops[l].begin();
			for(int i = 1; i <= p1; i++)
				iter_loop1++;
			Loop::iterator iter_loop2 = silLoops[l].begin();
			for(int i = 1; i <= p2; i++)
				iter_loop2++;

			LoopCirculator ccLoop1(&(silLoops[l]), iter_loop1);
			LoopCirculator ccLoop2(&(silLoops[l]), iter_loop2);
			
			while(ccLoop1 != ccLoop2){
				visShadows_iter->push_back(silPts[silLoops[l][*ccLoop1]]);
				ccLoop1++;
			};			
			visShadows_iter->push_back(silPts[silLoops[l][*ccLoop1]]);
*/
			if(p1 < p2){
				for(int i = p1; i <= p2; i++)
					visShadows_iter->push_back(silPts[silLoops[l][i]]);
			}
			else
			{
				for(int i = p1; i < silLoops[l].size(); i++)
					visShadows_iter->push_back(silPts[silLoops[l][i]]);
				for(int i = 0; i <= p2; i++)
					visShadows_iter->push_back(silPts[silLoops[l][i]]);			
			}
			//close polygon
			visShadows_iter->push_back(*(visShadows_iter->begin()));
/*			
			wxPoint* points = new wxPoint[ visShadows_iter->size()];
			int pIdx = 0;
			wxBrush wB;
			wB.SetColour(visShadows_iter->size(),250-visShadows_iter->size(),250);
			wxPen wP;
			wP.SetColour(visShadows_iter->size(),250-visShadows_iter->size(),250);
			wP.SetWidth(3);
			dc.SetBrush(wB);
			dc.SetPen(wP);
			for(visShadow_iter = visShadows_iter->begin(); 
				visShadow_iter != visShadows_iter->end(); ++visShadow_iter)
			{
					int xP = (int)CGAL::to_double((*visShadow_iter).x()); 
					int yP = (int)CGAL::to_double((*visShadow_iter).y());
					points[pIdx] = wxPoint(xP,yP);
					++pIdx;
			}
			dc.DrawPolygon( visShadows_iter->size(), points, 0, 0, wxWINDING_RULE);
	*/
		}
	}
	
	int numberOfLoops = visShadows.size();

	int* count = new int[numberOfLoops]; //array of loops sizes
	int* p_count = count;
	int n_points = 0; //total number of points
	for(visShadows_iter = visShadows.begin(); 
		visShadows_iter != visShadows.end(); ++visShadows_iter){
		*p_count = visShadows_iter->size();
		n_points += *p_count;
		p_count++;
	}
	
	wxPoint* points = new wxPoint[n_points];	
	int idx = 0; //index of the current point
	for(visShadows_iter = visShadows.begin(); 
		visShadows_iter != visShadows.end(); ++visShadows_iter){
		for(visShadow_iter = visShadows_iter->begin(); 
			visShadow_iter != visShadows_iter->end(); ++visShadow_iter){
			int xP = (int)CGAL::to_double((*visShadow_iter).x()); 
			int yP = (int)CGAL::to_double((*visShadow_iter).y());
			points[idx] = wxPoint(xP, yP);
			idx++;
		}
	}	
	
	wxBrush wB;
	wB.SetColour(255,200,150);
	wxPen wP;
	wP.SetColour(0,255,0);
	wP.SetWidth(3);
	dc.SetBrush(wB);
	dc.SetPen(wP);

	dc.DrawPolyPolygon( numberOfLoops, count,  points, 0, 0, wxWINDING_RULE);

	
//	myPM.insert (shadowCurves.begin(),shadowCurves.end());
/*	
	myPM.is_valid();
	std::cout << "no of faces: " << myPM.number_of_faces() << std::endl;
	int nFaces = myPM.number_of_faces();
	std::cout << "no of halfedges: " << myPM.number_of_halfedges() << std::endl;
	int nHedges = myPM.number_of_halfedges();
	std::cout << "no of vertices: " << myPM.number_of_vertices() << std::endl;
	int nVerts = myPM.number_of_vertices();

	//draw all faces

	int faceCount = 0;
	dc.SetPen(*wxGREEN_PEN);
	
	for (fit = myPM.faces_begin(); fit != myPM.faces_end(); ++fit) 
	{
		if (fit->is_unbounded()) continue;
		//compute the number of points
		int noFacePoints = 0;
		cc =  fit->outer_ccb();
		do{
			noFacePoints += (cc->curve()).points() - 1;
			++cc;
		} while(++cc != fit->outer_ccb());
		
		wxPoint* points = new wxPoint[ noFacePoints];
		int pIdx = 0;
		wxBrush wB;
		wxPen wP;
		wP.SetColour(faceCount,250-faceCount,250);
		faceCount += 60;
		//wB = *wxTRANSPARENT_BRUSH;
		wB.SetColour(faceCount,250-faceCount,250);
		dc.SetBrush(wB);
		dc.SetPen(wP);
		cc =  fit->outer_ccb();
		do{
			if(*((cc->curve()).begin()) == cc->target()->point())
				for(pit = (cc->curve()).begin(); pit != (cc->curve()).end(); ++pit)
				{
					int xP = (int)CGAL::to_double((*pit).x()); 
					int yP = (int)CGAL::to_double((*pit).y());
					points[pIdx] = wxPoint(xP,yP);
					++pIdx;
				}
			else
				for(pit = (cc->curve()).rbegin(); pit != (cc->curve()).rend(); ++pit)
				{
					int xP = (int)CGAL::to_double((*pit).x()); 
					int yP = (int)CGAL::to_double((*pit).y());
					points[pIdx] = wxPoint(xP,yP);
					++pIdx;
				}
			--pIdx;		
			++cc;
		} while(++cc != fit->outer_ccb());
		dc.DrawPolygon( noFacePoints, points, 0, 0, wxWINDING_RULE);
		delete[] points;
	}
 */
	//draw vertices
/*	
	Pmwx::Vertex_iterator vit;
	dc.SetPen(*wxRED_PEN);
	for(vit = myPM.vertices_begin(); vit != myPM.vertices_end(); vit++)
	{
		dc.DrawPoint((int)CGAL::to_double((vit->point()).x()), (int)CGAL::to_double((vit->point()).y()));
	}
*/	
}


void ShaderSVG::drawShade(wxSVGFileDC& dc, std::pair<Particles* const, int>& pSil_pair, std::pair<Particles* const, int>& pSha_pair)
{
	//get silhouette visibility
//	Particles* p_collectionSil = pSil_pair.first;
	
//	ParticleVisibility *silVis = dynamic_cast<ParticleVisibility *>(p_collectionSil->getAttributeGeneric("ParticleVisibility"));
//	if(!silVis) return;
	
	//get shadow visibility
//	Particles* p_collectionSha = pSha_pair.first;
	
//	ParticleVisibility *shaVis = dynamic_cast<ParticleVisibility *>(p_collectionSha->getAttributeGeneric("ParticleVisibility"));
//	if(!shaVis) return;
	
	//get camera position
//	gmVector3* cameraPosition = viewDep->getCameraPosition();
	
	computePlanarMap(dc, pSil_pair, pSha_pair);
	
}

bool ShaderSVG::findClosestParticle(const Point_2& sourcePt, 
									const std::vector<Loop>& loops, 
									std::pair<Particles* const, int>& pSil_pair, 
									Point_2& targetPt,
									int& p_index, //the closest particle index
									int& l_index) //the corresponding loop 
{
	//get silhouette visibility
	Particles* p_collection = pSil_pair.first;
	int found = pSil_pair.second;
	ShaderContour *sc = dynamic_cast<ShaderContour *>(p_collection->shaders[found]);
	ParticleVisibility *vis = dynamic_cast<ParticleVisibility *>(p_collection->getAttributeGeneric("ParticleVisibility"));
	if(!vis) return false;

	std::vector<Loop>::const_iterator loops_iter; //iterator over loops
	Loop::const_iterator loop_iter;	//iterator over a loop specified by loops_iter
	NT dist = gmGOOGOL;
	NT dist2;
	unsigned int j = 0;
	for(loops_iter = loops.begin(); loops_iter != loops.end(); ++loops_iter)
	{
		unsigned int l = 0;
		for(loop_iter = loops_iter->begin(); loop_iter != loops_iter->end(); ++loop_iter)
		{
			if(vis->getVisibility(*loop_iter) == 1){
				dist2 = CGAL::squared_distance(sourcePt,silPts[*loop_iter]);
				if(dist2 < dist){
					dist = dist2;
					targetPt = silPts[*loop_iter];
					p_index = l;
					l_index = j;
				}
			}
			l++;
		}
		j++;
	}
	return true;
}

bool ShaderSVG::findClosestPoint(const Point_2& sourcePt,  Point_2& targetPt)
{
	bool found = true;
	X_monotone_curve_2::const_iterator pit;
	Pmwx::Ccb_halfedge_circulator cc;
	Pmwx::Halfedge_handle hh;
	Pmwx::Locate_type lt;
	Pmwx::Edge_iterator eit;
	
	hh = myPM.locate(sourcePt,lt);
	NT dist2; 
	switch ( lt ){
		case 1: //VERTEX
			targetPt = Point_2((hh->target()->point()).x(), (hh->target()->point()).y(),-1);			
			break;
						
		case 2 : //EDGE:
			targetPt = *((hh->curve()).begin());
			dist2 = CGAL::squared_distance(sourcePt, targetPt);
			for(pit = (hh->curve()).begin(); pit != (hh->curve()).end(); ++pit){
				NT dist2next = CGAL::squared_distance(sourcePt, *pit);
				if(dist2next < dist2){
					dist2 = dist2next;
					targetPt = *pit;
				}
			}
			break;
						
		default : 			
			/*
			cc =  hh->face()->outer_ccb();
			targetPt = *((cc->curve()).begin());
			dist2 = CGAL::squared_distance(sourcePt, *((cc->curve()).begin()));
			do{
				for(pit = (cc->curve()).begin(); pit != (cc->curve()).end(); ++pit){
					NT dist2next = CGAL::squared_distance(sourcePt, *pit);
					if(dist2next < dist2){
						dist2 = dist2next;
						targetPt = *pit;
						}
					}
					++cc;
				} while(++cc !=  hh->face()->outer_ccb());
			 */
			targetPt = *(((myPM.edges_begin())->curve()).begin());
			dist2 = CGAL::squared_distance(sourcePt, targetPt);
			for(eit = myPM.edges_begin(); eit != myPM.edges_end(); eit++)
			{
				for(pit = (eit->curve()).begin(); pit != (eit->curve()).end(); ++pit){
					NT dist2next = CGAL::squared_distance(sourcePt, *pit);
					if(dist2next < dist2){
						dist2 = dist2next;
						targetPt = *pit;
					}
				}
			}
			break;
						
		//case 6 : //UNBOUNDED_FACE
		//default:
		//	found = false;
		//	break;
	} 
	return found;			
}

void ShaderSVG::computeVisualEvents(std::pair<Particles* const, int>& pSil_pair, std::pair<Particles* const, int>& pSha_pair)
{
	std::vector<Loop>::iterator loops_iter; //iterator over loops
												  //get silhouette visibility
	Particles* p_collection = pSil_pair.first;
	int found = pSil_pair.second;
	ShaderContour *sc = dynamic_cast<ShaderContour *>(p_collection->shaders[found]);
	int isVisble1 =1;
	int isVisble2 =1;
	ParticleVisibility *vis = dynamic_cast<ParticleVisibility *>(p_collection->getAttributeGeneric("ParticleVisibility"));
	ParticlePosition *pPos = dynamic_cast<ParticlePosition *>(p_collection->getAttributeGeneric("ParticlePosition"));
	ParticleNormal *pNormal = dynamic_cast<ParticleNormal *>(p_collection->getAttributeGeneric("ParticleNormal"));
	if(!vis) return;
	
	for(loops_iter = silLoops.begin(); loops_iter != silLoops.end(); ++loops_iter)
	{
		//initialize circulators
		LoopCirculator loopCc(&(*loops_iter));
		LoopCirculator loopCcNext(&(*loops_iter));
		LoopCirculator loopCcInit(&(*loops_iter));
		++loopCcNext;
		//find  self-intersections
		do{
			isVisble1 = vis->getVisibility(*loopCc);
			isVisble2 = vis->getVisibility(*loopCcNext);
			if(((isVisble1 == 0) && (isVisble2 == 1)) ||
			   ((isVisble1 == 1) && (isVisble2 == 0))){
				//find the intersection point
				gmVector3 p1 = pPos->getPosition(*loopCc);
				gmVector3 p2 = pPos->getPosition(*loopCcNext);
				gmVector3 cP = *(viewDep->getCameraPosition());
				gmVector3 n1 = pNormal->getNormal(*loopCc);
				gmVector3 n2 = pNormal->getNormal(*loopCcNext);
				double dot1 = dot(cP,n1);
				double dot2 = dot(cP,n2);
				if(dot1 - dot2 == 0) break;
				gmVector3 newPos = (1/(dot1 - dot2)) * (dot1 * p1 - dot2 * p2); 
				//silPts[*loopCc]
				//silPts[*loopCcNext]
			}
			++loopCc;
			++loopCcNext;
		}while(loopCc != loopCcInit);
	}
			
	
}

#endif

