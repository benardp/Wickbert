/**
 * Implementation of an SVG shader
 * @file SVGshader.cpp
 * @date 04/08/2006
 * @author Matei N. Stroila
 * @remarks
 */

#ifdef WB_USE_CGAL

//#define _DEBUG

#include "SVGshader.h"
#include "boost/multi_array.hpp"
#include "Particles/Shaders/ShaderSilhouetteContour.h"
#include "Particles/Behaviors/SilhouetteAdhesion.h"
#include "Particles/Shaders/ShaderSuggestiveContour.h"
#include "Particles/Behaviors/SuggestiveContourAdhesion.h"
#include "Particles/Shaders/ShaderShadowContour.h"
#include "Particles/Behaviors/ShadowAdhesion.h"
#include "Particles/Shaders/ShaderSpecularContour.h"
#include "Particles/Behaviors/SpecularAdhesion.h"
#include "Particles/Shaders/ShaderParabolicsContour.h"
#include "Particles/Behaviors/ParabolicsAdhesion.h"
#include "Particles/Shaders/ShaderDeCarloSuggestiveContour.h"
#include "Particles/Behaviors/DeCarloSuggestiveContourAdhesion.h"


#include "Particles/Attributes/ParticlePosition.h"
#include "Particles/Attributes/ParticleVisibility.h"
#include "Particles/Attributes/ImplicitInterrogator.h"
#include "Particles/Shaders/FindFirstIntersectionWithImplicitsShader.h"

#include <iterator>

typedef  CGAL::Const_circulator_from_container<Loop> LoopCirculator;

namespace
{
void convertOpenGLCoordinatesToNormalizedCoordinates(int viewport[4], double x, double y, double * resX, double * resY)
{
    *resX=(x-viewport[0])/(viewport[2]-viewport[0]);
    *resY=(y-viewport[1])/(viewport[3]-viewport[1]);
}



void convertNormalizedCoordinatesToOpenGLCoordinates(int viewport[4], double x, double y, double * resX, double * resY)
{
    *resX=(x*(viewport[2]-viewport[0]))+viewport[0];
    *resY=(y*(viewport[3]-viewport[1]))+viewport[1];
}

void convertNormalizedCoordinatesToDCCoordinates(const wxDC & dc, double x, double y, double * resX, double * resY)
{
    wxCoord width, height;
    dc.GetSize(&width, &height);
    *resX=(x*width);
    *resY=((1.0-y)*height);
}

void convertDCCoordinatesToNormalizedCoordinates(const wxDC & dc, double x, double y, double * resX, double * resY)
{
    wxCoord width, height;
    dc.GetSize(&width, &height);
    *resX=x/width;
    *resY=1.0-y/height;

}


void convertWorldPositionToOpenGLCoordinates(	int viewport[4], double projection[16] ,double modelview[16],
const gmVector3 & worldPosition, double * openGLX, double * openGLY, double * openGLZ)
{

    gluProject(
                worldPosition[0],
            worldPosition[1],
            worldPosition[2],
            modelview,
            projection,
            viewport,
            openGLX,
            openGLY,
            openGLZ
            );
}

void convertOpenGLCoordinatesToWorldPosition(int viewport[4], double projection[16], double modelview[16],
double x, double y, double depth, gmVector3 * point)
{
    gluUnProject( x,y,depth, modelview,
                  projection,
                  viewport,
                  &((*point)[0]),
            &((*point)[1]),
            &((*point)[2]) );
}


void insertNewImplicitFromInterrogators(ParticleStuff* particles, std::vector<Implicit*>* implicits)
{
    assert(implicits);
    if (!particles)
        return;
    std::vector<ImplicitInterrogator*> interrogators;
    particles->ps->getAllAttributesByType<ImplicitInterrogator>(&interrogators);
    for (unsigned int i=0;i<interrogators.size();++i)
    {
        Implicit* imp=(interrogators[i])->getImplicit();
        if (!imp)
            continue;
        if (find(implicits->begin(),implicits->end(),imp)==implicits->end())
            implicits->push_back(imp);
    }
}

template<class shader, class behavior>
void findShaderAndBehavior(const ParticleSystems * ps, const shader ** resShader, const behavior ** resBehavior)
{
    *resShader=0;
    *resBehavior=0;
    for (unsigned int s=0;s<ps->size();++s)
    {
        for(unsigned int i = 0; i < (*ps)[s]->particles.size(); i++)
        {
            Particles* currentPCollection = (*ps)[s]->particles[i];
            *resShader=currentPCollection->getShaderByType<shader>();

            if (*resShader)
            {
                *resBehavior=currentPCollection->getBehaviorByType<behavior>();
                if (*resBehavior)
                    return;
            }
        }
    }
}
}

wxBrush myTRANSPARENT_BRUSH; // needed to avoid const v. & reference with wxw2.8.0 -jch

SVGshader::SVGshader()
{
    myPs = 0;
    drawSilhouettesFlag = false;
    drawShadowsFlag = false;
    drawSpecularsFlag = false;
    drawSuggestivesFlag = false;
    drawDeCarloSuggestivesFlag = false;
    drawParabolicsFlag = false;

    myTRANSPARENT_BRUSH = *wxTRANSPARENT_BRUSH;
}

void SVGshader::setParticleSystems(ParticleSystems *ps)
{
    myPs = ps;
}


void SVGshader::drawSVG(wxDC& dc)
{
    if(!myPs) return;
    if(myPs->size() == 0 ) return;
    //initialize the canvas
    dc.SetBackground(*wxWHITE_BRUSH);
    dc.Clear();

    //look for all contours in the first particle system and project them
    //onto the view plane
    findContoursAndBehaviors();
    projectContoursOnNormalizedCoords();

    //draw silhouettes
    wxBrush wB;
    wxPen wP;
    if(shaderSilhouette){
        if(shaderSilhouette->IsShaded()){
            std::vector<short> brushColor;
            shaderSilhouette->getBrushColor(brushColor);
            wB.SetColour(brushColor[0],brushColor[1],brushColor[2]);
        }
        else
            wB = *wxTRANSPARENT_BRUSH;
        //set pen
        std::vector<short> penColor;
        shaderSilhouette->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderSilhouette->getPenWidth());
        drawSVGContours(dc, wP, wB, shaderSilhouette, silPts);
    }
    //draw shadows
    if(shaderShadow){
        //set pen
        std::vector<short> penColor;
        shaderShadow->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderShadow->getPenWidth());
        //draw contours
        drawSVGContours(dc, wP, myTRANSPARENT_BRUSH, shaderShadow, shadPts);
    }

    //draw suggestive contours
    if(shaderSuggestive){
        //set pen
        std::vector<short> penColor;
        shaderSuggestive->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderSuggestive->getPenWidth());
        //draw contours
        drawSVGContours(dc, wP, myTRANSPARENT_BRUSH, shaderSuggestive, suggPts);
    }

    //draw expensive suggestive contours
    if(shaderDeCarloSuggestive){
        //set pen
        std::vector<short> penColor;
        shaderDeCarloSuggestive->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderDeCarloSuggestive->getPenWidth());
        //draw contours
        drawSVGContours(dc, wP, myTRANSPARENT_BRUSH, shaderDeCarloSuggestive, deCarloSuggPts);
    }


    //draw parabolic contours
    if(shaderParabolic){
        //set pen
        std::vector<short> penColor;
        shaderParabolic->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderParabolic->getPenWidth());
        //draw contours
        drawSVGContours(dc, wP, myTRANSPARENT_BRUSH, shaderParabolic, parabPts);
    }
    //draw shaded shadows
    if(shaderSilhouette && shaderShadow && shaderShadow->IsShaded()){
        std::vector<short> shadeColor;
        shaderShadow->getShadeColor(shadeColor);
        wB.SetColour(shadeColor[0],shadeColor[1],shadeColor[2]);
        dc.SetBrush(wB);
        std::vector<short> penColor;
        shaderShadow->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderShadow->getPenWidth());
        dc.SetPen(wP);
        drawCycles(dc,shaderSilhouette,shaderShadow, shadPts );
    }

    //draw specular highlights
    if(shaderSilhouette && shaderSpecular && shaderSpecular->IsShaded()){
        std::vector<short> brushColor;
        shaderSpecular->getBrushColor(brushColor);
        wB.SetColour(brushColor[0],brushColor[1],brushColor[2]);
        dc.SetBrush(wB);
        //set pen
        std::vector<short> penColor;
        shaderSpecular->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderSpecular->getPenWidth());
        dc.SetPen(wP);
        drawCycles(dc,shaderSilhouette,shaderSpecular,specPts);
    }
}

void SVGshader::drawSVGContours(wxDC& dc, wxPen& pen, wxBrush& brush,  const ShaderContour* shaderContour, const std::vector<Point_2>& contPts, const bool drawContoursOnly)
{
    //get contour loops
    Particles* p_collection = shaderContour->ps;
    std::vector<Loop> loops = shaderContour->getLoops();
    //if no loops present, return
    if(loops.size() == 0) return;

    ParticlePosition *position = p_collection->getAttributeByType<ParticlePosition>();
    if(!position) return;

    std::vector<Loop>::const_iterator loops_iter; //iterator over loops
    Loop::const_iterator loop_iter, loop_iter_next;	//iterator over a loop specified by loops_iter

    if(!drawContoursOnly){
        ///////////////////
        //Draw shaded areas
        ///////////////////

        dc.SetBrush(brush);
        dc.SetPen(*wxTRANSPARENT_PEN);
        int* count = new int[loops.size()]; //array of loops sizes
        int* p_count = count;
        int n_points = 0; //total number of points
        for(loops_iter = loops.begin(); loops_iter != loops.end(); ++loops_iter){
            *p_count = loops_iter->size();
            n_points += *p_count;
            p_count++;
        }

        wxPoint* points = new wxPoint[n_points];
        int i = 0; //index of the current point
        for(loops_iter = loops.begin(); loops_iter != loops.end(); ++loops_iter){
            for(loop_iter = loops_iter->begin(); loop_iter != loops_iter->end(); ++loop_iter){
                double wx, wy;
                convertNormalizedCoordinatesToDCCoordinates(	dc,
                                                                CGAL::to_double(contPts[*loop_iter].x()),CGAL::to_double(contPts[*loop_iter].y()),
                        &wx,&wy);
                points[i] = wxPoint((int)wx, (int)wy);
                i++;
            }
        }
        if(loops.size() > 1)
            dc.DrawPolyPolygon( loops.size(), count,  points, 0, 0, wxWINDING_RULE);
        else
            dc.DrawPolygon( count[0],  points, 0, 0, wxWINDING_RULE);
        delete [] points;
        delete[] count;
    }

    ////////////////////
    //Draw contour loops
    ////////////////////
    ParticleVisibility *vis = p_collection->getAttributeByType<ParticleVisibility>();
    int isVisble1 =1;
    int isVisble2 =1;
    if(!vis) return;

    short invWidth = shaderContour->getInvPenWidth();
    wxPen wP;
    if(invWidth > 0){
        std::vector<short> invColor;
        shaderContour->getInvisibleContourColor(invColor);
        wP.SetColour(invColor[0], invColor[1], invColor[2]); //pen for invisible parts
        wP.SetWidth(invWidth);
    }
    else
        wP = *wxTRANSPARENT_PEN; //pen for invisible parts

    dc.SetBrush(*wxTRANSPARENT_BRUSH);
    dc.SetPen(pen);
    for(loops_iter = loops.begin(); loops_iter != loops.end(); ++loops_iter){
        loop_iter_next = loops_iter->begin();
        ++loop_iter_next;
        for(loop_iter = loops_iter->begin(); loop_iter != loops_iter->end(); ++loop_iter)
        {
            if(loop_iter_next == loops_iter->end()) loop_iter_next = loops_iter->begin();
            if(shaderContour->useVisibility){
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
            double x0, y0, x1, y1;
            convertNormalizedCoordinatesToDCCoordinates(	dc,
                                                            CGAL::to_double(contPts[*loop_iter].x()),CGAL::to_double(contPts[*loop_iter].y()),
                    &x0, &y0);
            convertNormalizedCoordinatesToDCCoordinates(	dc,
                                                            CGAL::to_double(contPts[*loop_iter_next].x()),CGAL::to_double(contPts[*loop_iter_next].y()),
                    &x1, &y1);

            dc.DrawLine((int)x0, (int)y0, (int)x1, (int)y1);
            ++loop_iter_next;
        }
    }
}

bool SVGshader::projectLoopsOnNormalizedCoords(const ShaderContour* sc, std::vector<Point_2>& pts)
{
    pts.clear();
    //get contour loops
    const std::vector<Loop>& loops = sc->getLoops();

    ParticlePosition *position = sc->ps->getAttributeByType<ParticlePosition>();
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


    //project the Silhouette loops
    int no_particles = sc->ps->size();
    for(int i=0;i<no_particles;i++)
    {
        double openGLX;
        double openGLY;
        double openGLZ;

        convertWorldPositionToOpenGLCoordinates(	viewport,projection, modelview, position->getPosition(i),
                                                    &openGLX,&openGLY, &openGLZ);
        double nX,nY;
        convertOpenGLCoordinatesToNormalizedCoordinates(viewport,openGLX, openGLY, &nX,&nY);

        Point_2 pt(nX, nY); //corresponding screen space position
        pts.push_back(pt);
    }
    return true;
}


void SVGshader::drawCycles(wxDC& dc, 
                           const ShaderContour* shaderSil,
                           const ShaderContour* shaderContour,
                           const std::vector<Point_2>& contPts)
{

    std::vector<VectorOfPoint_2> visContours;
    std::vector<EndsOfContour>  contourEnds;

    computeVisibleContours(shaderSil,
                           shaderContour,
                           contPts,
                           visContours,
                           contourEnds);

    //add visible silhouettes
    //get silhouette visibility
    ParticleVisibility* vis = shaderSil->ps->getAttributeByType<ParticleVisibility>();

    const std::vector<Loop>& silLoops = shaderSil->getLoops();

    std::vector<VectorOfPoint_2>::iterator visShadows_iter = visContours.begin();
    VectorOfPoint_2::iterator visShadow_iter;
    std::vector<EndsOfContour>::iterator ends_iter;
    for(ends_iter = contourEnds.begin(); ends_iter != contourEnds.end();
        ends_iter++, visShadows_iter++)
    {
        if(ends_iter->first.first == -1)
            continue;

        if(ends_iter->first.first == ends_iter->second.first) //same silhouette loop
        {
            int l = ends_iter->first.first;
            int p1 = ends_iter->second.second;
            int p2 = ends_iter->first.second;

            if(p1 < p2){
                for(int i = p1; i <= p2; i++)
                    if(vis->getVisibility(silLoops[l][i]) == 1)
                        visShadows_iter->push_back(silPts[silLoops[l][i]]);
            }
            else
            {
                for(unsigned int i = p1; i < silLoops[l].size(); i++)
                    if(vis->getVisibility(silLoops[l][i]) == 1)
                        visShadows_iter->push_back(silPts[silLoops[l][i]]);
                for(int i = 0; i <= p2; i++)
                    if(vis->getVisibility(silLoops[l][i]) == 1)
                        visShadows_iter->push_back(silPts[silLoops[l][i]]);
            }
            //close polygon
            //visShadows_iter->push_back(*(visShadows_iter->begin()));
        }
    }

    int numberOfLoops = visContours.size();

    int* count = new int[numberOfLoops]; //array of loops sizes
    int* p_count = count;
    int n_points = 0; //total number of points
    for(visShadows_iter = visContours.begin();
        visShadows_iter != visContours.end(); ++visShadows_iter){
        *p_count = visShadows_iter->size();
        n_points += *p_count;
        p_count++;
    }

    wxPoint* points = new wxPoint[n_points];
    int idx = 0; //index of the current point
    for(visShadows_iter = visContours.begin();
        visShadows_iter != visContours.end(); ++visShadows_iter){
        for(visShadow_iter = visShadows_iter->begin();
            visShadow_iter != visShadows_iter->end(); ++visShadow_iter){
            double xP, yP;
            convertNormalizedCoordinatesToDCCoordinates(	dc,
                                                            CGAL::to_double((*visShadow_iter).x()), CGAL::to_double((*visShadow_iter).y()),
                                                            &xP,&yP);
            points[idx] = wxPoint((int)xP, (int)yP);
            idx++;
        }
    }

    if(numberOfLoops > 1)
        dc.DrawPolyPolygon( numberOfLoops, count,  points, 0, 0, wxWINDING_RULE);
    else
        dc.DrawPolygon( count[0],  points, 0, 0, wxWINDING_RULE);

}

void SVGshader::computeVisibleContours(const ShaderContour* shaderSil, 
                                       const ShaderContour* shaderContour,
                                       const std::vector<Point_2> contour2DPts,
                                       std::vector<VectorOfPoint_2>& visContours,
                                       std::vector<EndsOfContour>& contourEnds)
{	
    //new visibility chaining:
    //for silhouettes we don't need to find the closest silhouette particles
    bool sameContours = false;
    if(shaderContour->ps == shaderSil->ps) sameContours = true;


    ParticleVisibility *vis = shaderContour->ps->getAttributeByType<ParticleVisibility>();
    if(!vis) return;


    std::vector<Loop>::const_iterator loops_iter; //iterator over loops
    const std::vector<Loop>& contourLoops = shaderContour->getLoops();
    for(loops_iter = contourLoops.begin(); loops_iter != contourLoops.end(); loops_iter++)
    {
        //we only want to keep really connected stuff.
        //so at least 3 particles should be connected in a loop.
        if (loops_iter->size()<3)
            continue;

        //initialize circulators
        LoopCirculator loopCc(&(*loops_iter));
        LoopCirculator loopCcNext(&(*loops_iter));
        LoopCirculator loopCcInit(&(*loops_iter));
        ++loopCcNext;

        //we go to first invisible->visible event, if we end up in a circle, we stop...
        bool allHasSameVisibility=false;
        while(!((vis->getVisibility(*loopCc) == 0)&&(vis->getVisibility(*loopCcNext)==1)))
        {
            ++loopCc; ++loopCcNext;
            if (loopCc==loopCcInit)
            {
                allHasSameVisibility=true;
                break;
            }
        }
        if (allHasSameVisibility)
        {
            if (vis->getVisibility(*loopCcInit) == 0)
                continue;
        }

        loopCcInit=loopCc;



        VectorOfPoint_2 loop_pts;
        int p_idx, l_idx;

        EndsOfContour s_ends;
        ContourEnd s1_end, s2_end;
        std::vector<EndsOfContour>::iterator ends_iter;
        bool onEnd = false;
        Point_2 closestPoint;
        int isVisble1 =0;
        int isVisble2 =1;

        do{
            //this should not be necessary here...
            isVisble1 = vis->getVisibility(*loopCc);
            isVisble2 = vis->getVisibility(*loopCcNext);


            onEnd = false;
            //visible event: invisible -> visible
            if((isVisble1 == 0) && (isVisble2 == 1)){
                //find closest silhouette particle
                if(!sameContours){
                    findClosestParticle(contour2DPts[*loopCcNext], shaderSil, closestPoint,  p_idx, l_idx);
                    loop_pts.push_back(	closestPoint );
                    std::cout << "closest sil particle on opening part:" << closestPoint << std::endl;
                    s1_end = std::make_pair(l_idx, p_idx);
                }
                else
                {
                    loop_pts.push_back(contour2DPts[*loopCc]);
                }
                ++loopCc;
                ++loopCcNext;
            }
            //visible event: visible -> invisible
            if((isVisble1 == 1) && (isVisble2 == 0)){
                loop_pts.push_back(contour2DPts[*loopCc]);
                //std::cout << "closing part:" << contour2DPts[*loopCc] << std::endl;
                //find closest silhouette particle
                if(!sameContours){
                    findClosestParticle(contour2DPts[*loopCc], shaderSil, closestPoint, p_idx, l_idx);
                    loop_pts.push_back(	closestPoint );
                    //std::cout << "closest sil particle:" << closestPoint << std::endl;
                    s2_end = std::make_pair(l_idx, p_idx);
                }
                else
                {
                    loop_pts.push_back(contour2DPts[*loopCcNext]);
                }
                ++loopCc;
                ++loopCcNext;
                onEnd = true;
            }
            //on invisible part
            while((isVisble1 == 0) && (isVisble2 == 0)){
                ++loopCc;
                ++loopCcNext;
                isVisble1 = vis->getVisibility(*loopCc);
                isVisble2 = vis->getVisibility(*loopCcNext);
                if(loopCc == loopCcInit) break;
            }
            //on visible part
            while((isVisble1 == 1) && (isVisble2 == 1)){
                loop_pts.push_back(contour2DPts[*loopCc]);
                ++loopCc;
                ++loopCcNext;
                isVisble1 = vis->getVisibility(*loopCc);
                isVisble2 = vis->getVisibility(*loopCcNext);
                if(loopCc == loopCcInit) {
                    loop_pts.push_back(contour2DPts[*loopCc]);
                    s1_end = std::make_pair(-1, -1);
                    s2_end = std::make_pair(-1, -1);
                    onEnd = true;
                    break;
                }
            }
            if(loop_pts.size()>1 && onEnd)
            {
                visContours.push_back(loop_pts);
                s_ends = std::make_pair(s1_end, s2_end);
                contourEnds.push_back(s_ends);
            }
            if(onEnd) loop_pts.clear();
        }while(loopCc != loopCcInit);
    }



    //for silhouettes we don't need to find the closest silhouette particles
    //bool sameContours = false;
    //if(shaderContour->ps == shaderSil->ps) sameContours = true;
    //
    //int isVisble1 =1;
    //int isVisble2 =1;
    //ParticleVisibility *vis = shaderContour->ps->getAttributeByType<ParticleVisibility>();
    //if(!vis) return;
    //VectorOfPoint_2 loop_pts;
    //std::vector<Loop>::const_iterator loops_iter; //iterator over loops
    //int p_idx, l_idx;
    //
    //EndsOfContour s_ends;
    //ContourEnd s1_end, s2_end;
    //std::vector<EndsOfContour>::iterator ends_iter;
    //bool onEnd = false;
    //Point_2 closestPoint;
    //
    //const std::vector<Loop>& contourLoops = shaderContour->getLoops();
    //for(loops_iter = contourLoops.begin(); loops_iter != contourLoops.end(); loops_iter++)
    //{
    //	//initialize circulators
    //	LoopCirculator loopCc(&(*loops_iter));
    //	LoopCirculator loopCcNext(&(*loops_iter));
    //	LoopCirculator loopCcInit(&(*loops_iter));
    //	++loopCcNext;
    //	//move to the first visibility event
    //	while((vis->getVisibility(*loopCcInit) == vis->getVisibility(*loopCcNext))
    //		&& (loopCcNext != loopCcInit )){
    //		++loopCcNext;
    //		++loopCc;
    //	}
    //	if(loopCcNext == loopCcInit && vis->getVisibility(*loopCcInit) == 0) continue; //no visible parts, move to the next loop
    //	//if on invisible part move to the next visible part
    //	if(vis->getVisibility(*loopCcNext) == 0){
    //		loopCcInit = loopCcNext;
    //		++loopCcNext;
    //		++loopCc;
    //		while((vis->getVisibility(*loopCcInit) == vis->getVisibility(*loopCcNext)) && (loopCcNext != loopCcInit )){
    //			++loopCcNext;
    //			++loopCc;
    //		}
    //	}
    //	//got to the first visible event: invisible -> visible
    //	loopCcInit = loopCc;
    //	do{
    //		onEnd = false;
    //		isVisble1 = vis->getVisibility(*loopCc);
    //		isVisble2 = vis->getVisibility(*loopCcNext);
    //		//visible event: invisible -> visible
    //		if((isVisble1 == 0) && (isVisble2 == 1)){
    //			//find closest silhouette particle
    //			if(!sameContours){
    //				findClosestParticle(contour2DPts[*loopCcNext], shaderSil, closestPoint,  p_idx, l_idx);
    //				loop_pts.push_back(	closestPoint );
    //				std::cout << "closest sil particle on opening part:" << closestPoint << std::endl;
    //				s1_end = std::make_pair(l_idx, p_idx);
    //			}
    //			else
    //			{
    //				loop_pts.push_back(contour2DPts[*loopCc]);
    //			}
    //			++loopCc;
    //			++loopCcNext;
    //		}
    //		//visible event: visible -> invisible
    //		if((isVisble1 == 1) && (isVisble2 == 0)){
    //			loop_pts.push_back(contour2DPts[*loopCc]);
    //			//std::cout << "closing part:" << contour2DPts[*loopCc] << std::endl;
    //			//find closest silhouette particle
    //			if(!sameContours){
    //				findClosestParticle(contour2DPts[*loopCc], shaderSil, closestPoint, p_idx, l_idx);
    //				loop_pts.push_back(	closestPoint );
    //				//std::cout << "closest sil particle:" << closestPoint << std::endl;
    //				s2_end = std::make_pair(l_idx, p_idx);
    //			}
    //			else
    //			{
    //				loop_pts.push_back(contour2DPts[*loopCcNext]);
    //			}
    //			++loopCc;
    //			++loopCcNext;
    //			onEnd = true;
    //		}
    //		//on invisible part
    //		while((isVisble1 == 0) && (isVisble2 == 0)){
    //			++loopCc;
    //			++loopCcNext;
    //			isVisble1 = vis->getVisibility(*loopCc);
    //			isVisble2 = vis->getVisibility(*loopCcNext);
    //			if(loopCc == loopCcInit) break;
    //		}
    //		//on visible part
    //		while((isVisble1 == 1) && (isVisble2 == 1)){
    //			loop_pts.push_back(contour2DPts[*loopCc]);
    //			++loopCc;
    //			++loopCcNext;
    //			isVisble1 = vis->getVisibility(*loopCc);
    //			isVisble2 = vis->getVisibility(*loopCcNext);
    //			if(loopCc == loopCcInit) {
    //				loop_pts.push_back(contour2DPts[*loopCc]);
    //				s1_end = std::make_pair(-1, -1);
    //				s2_end = std::make_pair(-1, -1);
    //				onEnd = true;
    //				break;
    //			}
    //		}
    //		if(loop_pts.size()>1 && onEnd)
    //		{
    //			visContours.push_back(loop_pts);
    //			s_ends = std::make_pair(s1_end, s2_end);
    //			contourEnds.push_back(s_ends);
    //		}
    //		if(onEnd) loop_pts.clear();
    //	}while(loopCc != loopCcInit);
    //}


}

bool SVGshader::findClosestParticle(const Point_2& sourcePt, 
                                    const ShaderContour* shaderSil,
                                    Point_2& targetPt,
                                    int& p_index, //the closest particle index
                                    int& l_index) //the corresponding loop
{
    //get silhouette visibility
    ParticleVisibility *vis = shaderSil->ps->getAttributeByType<ParticleVisibility>();
    if(!vis) return false;

    std::vector<Loop>::const_iterator loops_iter; //iterator over loops
    Loop::const_iterator loop_iter;	//iterator over a loop specified by loops_iter
    double dist = gmGOOGOL;
    double dist2;
    unsigned int j = 0;

    const std::vector<Loop>& loops = shaderSil->getLoops();
    for(loops_iter = loops.begin(); loops_iter != loops.end(); ++loops_iter)
    {
        unsigned int l = 0;
        for(loop_iter = loops_iter->begin(); loop_iter != loops_iter->end(); ++loop_iter)
        {
            if(vis->getVisibility(*loop_iter) == 1){
                dist2 = CGAL::to_double(CGAL::squared_distance(sourcePt,silPts[*loop_iter]));
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


void SVGshader::builAndDrawPlanarMap(wxDC& dc)
{
    if(!myPs) return;
    if(myPs->size() == 0 ) return;
    //clear the planar map
    arr.clear();
    cdt.clear();

    //look for all contours in the first particle system and project them
    //onto the view plane
    findContoursAndBehaviors();
    projectContoursOnNormalizedCoords();

    std::vector<VectorOfPoint_2>::iterator visContours_iter;
    VectorOfPoint_2::iterator iterPoint1, iterPoint2;
    //insert shadows in planar map
    if(shaderSilhouette && shaderShadow){
        std::vector<VectorOfPoint_2>  visContour;
        std::vector<EndsOfContour>  contourEnds;
        computeVisibleContours(shaderSilhouette,
                               shaderShadow,
                               shadPts,
                               visContour,
                               contourEnds);
        std::cout << "no of vis contours:" << visContour.size() << std::endl;

        for(visContours_iter = visContour.begin(); visContours_iter != visContour.end(); ++visContours_iter){
            //data 1 means shadows for now
            CGAL::insert_curve (arr, Curve_2(Polyline_2(visContours_iter->begin(), visContours_iter->end()), 1));
            //insert constraints in the constrained triangulation
            iterPoint2 = visContours_iter->begin();
            iterPoint2++;
            for(iterPoint1 = visContours_iter->begin(); iterPoint2 != visContours_iter->end(); iterPoint1++, iterPoint2++)
            {
                cdt.insert_constraint(*iterPoint1, *iterPoint2);
            }
        }
    }
    //insert speculars in planar map
    if(shaderSilhouette && shaderSpecular){
        std::vector<VectorOfPoint_2>  visContour;
        std::vector<EndsOfContour>  contourEnds;
        computeVisibleContours(shaderSilhouette,
                               shaderSpecular,
                               specPts,
                               visContour,
                               contourEnds);
        //std::cout << "no of vis contours:" << visContour.size() << std::endl;

        for(visContours_iter = visContour.begin(); visContours_iter != visContour.end(); ++visContours_iter){
            //data 2 means specularities for now
            CGAL::insert_curve (arr, Curve_2(Polyline_2(visContours_iter->begin(), visContours_iter->end()), 2));
            //insert constraints in the constrained triangulation
            iterPoint2 = visContours_iter->begin();
            iterPoint2++;
            for(iterPoint1 = visContours_iter->begin(); iterPoint2 != visContours_iter->end(); iterPoint1++, iterPoint2++)
            {
                cdt.insert_constraint(*iterPoint1, *iterPoint2);
            }
        }
    }
    //insert silhouettes in planar map
    if(shaderSilhouette){
        std::vector<VectorOfPoint_2>  visContour;
        std::vector<EndsOfContour>  contourEnds;
        computeVisibleContours(shaderSilhouette,
                               shaderSilhouette,
                               silPts,
                               visContour,
                               contourEnds);
#ifdef _DEBUG
        std::cout << "no of vis contours:" << visContour.size() << std::endl;
#endif
        for(visContours_iter = visContour.begin(); visContours_iter != visContour.end(); ++visContours_iter){
            //data 0 means silhouettes for now
            CGAL::insert_curve (arr, Curve_2(Polyline_2(visContours_iter->begin(), visContours_iter->end()), 0));
            //insert constraints in the constrained triangulation
            iterPoint2 = visContours_iter->begin();
            iterPoint2++;
            for(iterPoint1 = visContours_iter->begin(); iterPoint2 != visContours_iter->end(); iterPoint1++, iterPoint2++)
            {
                cdt.insert_constraint(*iterPoint1, *iterPoint2);
            }
        }
    }
    assert(cdt.is_valid());
    assert(arr.is_valid());
    //draw planar map
    draw_arrangement(dc);

}

//draw constrained Delaunay triangulation
void SVGshader::buildAndDrawConstrainedDelaunay(wxDC& dc)
{
    if(cdt.number_of_faces() == 0) return;
    if(arr.number_of_faces() == 0) return;

    std::vector<Point_2> barycenters;

    dc.SetBackground(*wxWHITE_BRUSH);
    dc.SetPen(*wxBLACK_PEN);
    dc.SetBrush(*wxGREEN_BRUSH);
    CDT::Finite_faces_iterator ffiter;
    for( ffiter = cdt.finite_faces_begin();  ffiter != cdt.finite_faces_end(); ffiter++){

        std::vector<wxPoint> points;
        std::vector<Point_2> facePoints;


        double x, y;
        convertNormalizedCoordinatesToDCCoordinates(	dc,
                                                        CGAL::to_double((ffiter->vertex(0)->point()).x()), CGAL::to_double((ffiter->vertex(0)->point()).y()),
                                                        &x, &y);
        points.push_back(wxPoint((int) x, (int)y));
        facePoints.push_back(ffiter->vertex(0)->point());

        convertNormalizedCoordinatesToDCCoordinates(	dc,
                                                        CGAL::to_double((ffiter->vertex(1)->point()).x()), CGAL::to_double((ffiter->vertex(1)->point()).y()),
                                                        &x, &y);

        points.push_back(wxPoint((int) x, (int) y));
        facePoints.push_back(ffiter->vertex(1)->point());


        convertNormalizedCoordinatesToDCCoordinates(	dc,
                                                        CGAL::to_double((ffiter->vertex(2)->point()).x()), CGAL::to_double((ffiter->vertex(2)->point()).y()),
                                                        &x, &y);
        points.push_back(wxPoint((int) x, (int) y));
        facePoints.push_back(ffiter->vertex(2)->point());

        dc.DrawPolygon(3, &(points[0]));

        Kernel::Vector_2 v1 = facePoints[0] - CGAL::ORIGIN;
        Kernel::Vector_2 v2 = facePoints[1] - CGAL::ORIGIN;
        Kernel::Vector_2 v3 = facePoints[2] - CGAL::ORIGIN;
        Point_2 centroid = CGAL::ORIGIN + (v1 + v2 + v3)*(1.0/3.0);
        barycenters.push_back(centroid);
    }

    //compute classifyFaceRaySets and interiorPoints vectors
    GLdouble projection[16];
    glGetDoublev(GL_PROJECTION_MATRIX, projection);// Retrieve The Projection Matrix
    GLdouble modelview[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);// Retrieve The Modelview Matrix
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);


    //initialize classifyFaceRaySets
    classifyFaceRaySets.clear();
    std::vector< std::pair<gmVector3, gmVector3> > faceRays;
    for (Arrangement_2::Face_const_iterator fit = arr.faces_begin(); fit != arr.faces_end(); ++fit){
        classifyFaceRaySets.push_back(faceRays);
    }

    // Perform the point-location query.
    PointLocation pl(arr);
    Arrangement_2::Face_const_handle  f;
    //loop over the centroids
    static const int skip = 3;
    for(unsigned int  i=0; i < barycenters.size(); i++){
        //locate the ith centroid
        CGAL::Object obj = pl.locate (barycenters[i]);
        if (CGAL::assign (f, obj)) {
            // q is located inside a face
            if (f->is_unbounded()){
#ifdef _DEBUG
                std::cout << "point located inside the unbounded face." << std::endl;
#endif
            }
            else
            {
                int j =0;
                for (Arrangement_2::Face_const_iterator fit = arr.faces_begin(); fit != arr.faces_end(); ++fit, ++j){
                    if(f == fit){
                        std::pair<gmVector3, gmVector3> ray;
                        //computeFaceRay(ray,dc, CGAL::to_double( barycenters[i].x()), CGAL::to_double( barycenters[i].y()),
                        //	  projection, modelview, viewport);
                        double openGLx,openGLy;
                        convertNormalizedCoordinatesToOpenGLCoordinates(	viewport,
                                                                            CGAL::to_double(barycenters[i].x()), CGAL::to_double(barycenters[i].y()),
                                                                            &openGLx,&openGLy);

                        convertOpenGLCoordinatesToWorldPosition(viewport,projection,modelview,openGLx,openGLy,0,&(ray.first));
                        convertOpenGLCoordinatesToWorldPosition(viewport,projection,modelview,openGLx,openGLy,1,&(ray.second));

                        classifyFaceRaySets[j].push_back(ray);
                        break;
                    }
                }
            }
        }
        else
        {
#ifdef _DEBUG
            std::cout << "point located on edge or vertex." << std::endl;
#endif
        }
        i += skip;
    }//end loop over barycenters

}

//look for all contours in the first particle system and project them
//onto the view plane
void SVGshader::findContoursAndBehaviors()
{
    if(myPs->size() == 0 ) return;

    //find contour shaders and associated behavior in particle systems
    //anyway this is a big hack so let's do it at least right ;)
    findShaderAndBehavior<ShaderSilhouetteContour, SilhouetteAdhesion>(myPs, &shaderSilhouette,&behaviorSilhouette);
    findShaderAndBehavior<ShaderShadowContour, ShadowAdhesion>(myPs, &shaderShadow,&behaviorShadow);
    findShaderAndBehavior<ShaderSpecularContour, SpecularAdhesion>(myPs, &shaderSpecular,&behaviorSpecular);
    findShaderAndBehavior<ShaderSuggestiveContour, SuggestiveContourAdhesion>(myPs, &shaderSuggestive,&behaviorSuggestive);
    findShaderAndBehavior<ShaderParabolicsContour, ParabolicsAdhesion>(myPs, &shaderParabolic,&behaviorParabolic);
    findShaderAndBehavior<ShaderDeCarloSuggestiveContour, DeCarloSuggestiveContourAdhesion>(myPs, &shaderDeCarloSuggestive,&behaviorSuggestiveExp);
}
void SVGshader::projectContoursOnNormalizedCoords()
{
    if (shaderSilhouette)
        projectLoopsOnNormalizedCoords( shaderSilhouette, silPts);

    if (shaderShadow)
        projectLoopsOnNormalizedCoords( shaderShadow, shadPts);

    if (shaderSpecular)
        projectLoopsOnNormalizedCoords( shaderSpecular, specPts);

    if (shaderSuggestive)
        projectLoopsOnNormalizedCoords( shaderSuggestive, suggPts);

    if (shaderDeCarloSuggestive)
        projectLoopsOnNormalizedCoords( shaderDeCarloSuggestive, deCarloSuggPts);

    if (shaderParabolic)
        projectLoopsOnNormalizedCoords( shaderParabolic, parabPts);
}


void  SVGshader::findImplicits(std::vector<Implicit*>* implicits)
{
    //to be honest... I would have taken only the implicits coming from the involved particles
    //and this is what I am gonna do...
    if(myPs->size() == 0 ) return;

    //Surfaces * surfs=(*myPs)[0]->surfaces;
    //if (!surfs) return;

    //Surfaces::iterator iter = surfs->begin();
    //const Surfaces::iterator end= surfs->end();

    //for (;iter!=end;++iter)
    //{
    //	Implicit* imp=dynamic_cast<Implicit*>(*iter);
    //	if (imp)
    //		implicits->push_back(imp);
    //}

    insertNewImplicitFromInterrogators((ParticleStuff*)shaderSilhouette, implicits);
    insertNewImplicitFromInterrogators((ParticleStuff*)shaderShadow, implicits);
    insertNewImplicitFromInterrogators((ParticleStuff*)shaderSpecular, implicits);
    insertNewImplicitFromInterrogators((ParticleStuff*)shaderParabolic, implicits);
    insertNewImplicitFromInterrogators((ParticleStuff*)shaderSuggestive, implicits);

    if (implicits->size()>1)
        std::cerr<<"Error there are several implicit surfaces involved in this example... check the implicit interrogators\n";
}

//Draw the planar map 
//-----------------------------------------------------------------------------
// Draw all edges (curves) along a connected component
// boundary.
//
unsigned int SVGshader::addFaceVerticesInDCCoordinates(	const wxDC & dc, 
                                                        Arrangement_2::Ccb_halfedge_const_circulator circ,
                                                        std::vector<wxPoint>& points)
{
    Arrangement_2::Ccb_halfedge_const_circulator  curr = circ;
    Arrangement_2::Halfedge_const_handle          he;
    unsigned int size = 0;
    do
    {
        he = curr;
        if((he->source()->degree() == 1) || (he->target()->degree() == 1 )){
            ++curr;
            continue;
        }
        if(he->source()->point() == *(he->curve().begin())){
            Curve_2::const_iterator         pit;
            for(pit = he->curve().begin(); pit != he->curve().end(); pit++){
                double x,y;
                convertNormalizedCoordinatesToDCCoordinates(dc, CGAL::to_double( pit->x()), CGAL::to_double( pit->y()),&x,&y);
                points.push_back(wxPoint((int)x, (int)y));
                size++;
            }
        }
        else
        {
            Curve_2::const_reverse_iterator         pit;
            for(pit = he->curve().rbegin(); pit != he->curve().rend(); pit++){
                double x,y;
                convertNormalizedCoordinatesToDCCoordinates(dc, CGAL::to_double( pit->x()), CGAL::to_double( pit->y()),&x,&y);
                points.push_back(wxPoint((int)x, (int)y));
                size++;
            }
        }
        ++curr;
    } while (curr != circ);
    return size;
}

//-----------------------------------------------------------------------------
// Draw face.
//
void SVGshader::draw_face (Arrangement_2::Face_const_handle f, wxDC& dc, 
                           wxBrush& wB,
                           const bool randomColour)
{
    // Print the outer boundary.
    if (f->is_unbounded())
    {
#ifdef _DEBUG
        std::cout << "Unbounded face. " << std::endl;
#endif
        return;
    }

    std::vector<int> count;
    std::vector<wxPoint> points;
    unsigned int nbNewPoints;
    nbNewPoints=addFaceVerticesInDCCoordinates(dc, f->outer_ccb(), points);
    count.push_back(nbNewPoints);

    Arrangement_2::Hole_const_iterator  hole;
    for (hole = f->holes_begin(); hole != f->holes_end(); ++hole)
    {
        nbNewPoints=addFaceVerticesInDCCoordinates(dc, *hole, points);
        count.push_back(nbNewPoints);
    }


    ///////////////////
    //Draw shaded areas
    ///////////////////
    if(randomColour){
        short color[3] = {short(rand()%255),short(rand()%255), short(rand()%255)};
        wB.SetColour(color[0],color[1],color[2]);
        wxPen pen;
        pen.SetColour(color[0]+10,color[1]+10,color[2]+10);
        dc.SetPen(pen);
    }
    dc.SetBrush(wB);
    if(count.size() > 1)
        dc.DrawPolyPolygon( count.size(), &(count[0]),  &(points[0]), 0, 0, wxWINDING_RULE);
    else
        dc.DrawPolygon( count[0],  &(points[0]), 0, 0, wxWINDING_RULE);

}

//-----------------------------------------------------------------------------
// Print the given arrangement.
//
void SVGshader::draw_arrangement (wxDC& dc)
{
    //CGAL::CGAL_precondition (arr.is_valid());
    // Draw the arrangement faces.
    Arrangement_2::Face_const_iterator    fit;
#ifdef _DEBUG	
    std::cout << arr.number_of_faces() << " faces:" << std::endl;
#endif
    wxBrush wB;
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
        draw_face(fit, dc, wB, true);

    // Draw the arrangement vertices.
    Arrangement_2::Vertex_const_iterator  vit;
#ifdef _DEBUG	
    std::cout << arr.number_of_vertices() << " vertices:" << std::endl;
#endif	
    for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
    {
        if (vit->is_isolated()){
#ifdef _DEBUG
            std::cout << " - Isolated." << std::endl;
#endif
            dc.SetPen(*wxRED_PEN);
        }
        else{
            dc.SetPen(*wxGREEN_PEN);
#ifdef _DEBUG
            std::cout << " - degree " << vit->degree() << std::endl;
#endif
        }
        double x,y;

        convertNormalizedCoordinatesToDCCoordinates(	dc,
                                                        CGAL::to_double( vit->point().x()),
                                                        CGAL::to_double( vit->point().y()),
                                                        &x, &y);
        dc.DrawCircle( (int)x,  (int)y,  4 );
    }
#ifdef _DEBUG		
    std::cout << arr.number_of_edges() << " edges:" << std::endl;
#endif

}

void SVGshader::buildSkeleton(wxDC& dc)
{
    if(!myPs) return;
    if(myPs->size() == 0 ) return;
    if(arr.number_of_faces() == 0) return;

    //Build the skeleton of each arrangement faces.
    classifyFaceRaySets.clear();
    interiorPoints.clear();
    Arrangement_2::Face_const_iterator    fit;
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
        constructSkeleton(fit, dc);

}

//Construct skeleton of a face
void SVGshader::constructSkeleton (Arrangement_2::Face_const_handle f, wxDC& dc)
{
    std::vector<std::pair<gmVector3,gmVector3> > faceTestRays;
    std::vector<Kernel2::Point_2> interiorPointsOfFace;
    if (f->is_unbounded())
    {
        //add the empty vector of test points
        interiorPoints.push_back(interiorPointsOfFace);
        classifyFaceRaySets.push_back(faceTestRays);
#ifdef _DEBUG
        std::cout << "Unbounded face. " << std::endl;
#endif
        return;
    }
    // Instantiate the skeleton builder
    SsBuilder ssb ;

    ccbForSkeleton(f->outer_ccb(), ssb);

    Arrangement_2::Hole_const_iterator  hole;
    for (hole = f->holes_begin(); hole != f->holes_end(); ++hole)
    {
        ccbForSkeleton (*hole,ssb);
    }

    // Construct the skeleton
    boost::shared_ptr<Ss> ss;
    try
    {
        ss = ssb.construct_skeleton();
    }
    catch( std::exception const& e )
    {
        std::cout << "EXCEPTION THROWN (" << e.what() << ") during straight skeleton construction.";
        return;
    }
    if(!ss){
        std::cout << "Unable to construct skeleton." << std::endl;
        return;
    }
    
    GLdouble projection[16];
    glGetDoublev(GL_PROJECTION_MATRIX, projection);// Retrieve The Projection Matrix
    GLdouble modelview[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);// Retrieve The Modelview Matrix
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    Ss::Vertex_const_iterator  vit;
    for (vit = ss->vertices_begin(); vit != ss->vertices_end(); ++vit)
    {
        if(vit->is_skeleton()){
            dc.SetPen(*wxBLACK_PEN);

            std::pair<gmVector3, gmVector3> ray;

            double openGLx,openGLy;
            convertNormalizedCoordinatesToOpenGLCoordinates(	viewport,
                                                                CGAL::to_double( vit->point().x()), CGAL::to_double( vit->point().y()),
                                                                &openGLx,&openGLy);


            convertOpenGLCoordinatesToWorldPosition(viewport,projection,modelview,openGLx,openGLy,0,&(ray.first));
            convertOpenGLCoordinatesToWorldPosition(viewport,projection,modelview,openGLx,openGLy,1,&(ray.second));

            faceTestRays.push_back(ray);
            interiorPointsOfFace.push_back(vit->point());
            double wx,wy;
            convertNormalizedCoordinatesToDCCoordinates(	dc,
                                                            CGAL::to_double(vit->point().x()), CGAL::to_double(vit->point().y()),&wx,&wy);

            dc.DrawCircle( (int) wx,  (int) wy,  2 );
        }
    }
    //add the face test points
    classifyFaceRaySets.push_back(faceTestRays);
    interiorPoints.push_back(interiorPointsOfFace);
}

//void SVGshader::computeFaceRay(std::pair<gmVector3, gmVector3>& ray, wxDC& dc, const double vx, const double vy, 
//							     double* projection,   double* modelview,   int* viewport )
//{
//
//	double normalizedX, normalizedY;
//	convertDCCoordinatesToNormalizedCoordinates( dc,vx,vy,
//												&normalizedX, &normalizedY);
//	double x,y; 
//	convertNormalizedCoordinatesToOpenGLCoordinates(viewport,normalizedX, normalizedY,&x,&y);			
//	gmVector3 point;
//	convertOpenGLCoordinatesToWorldPosition(viewport, projection, modelview,x,y,1,&point);
//	gmVector3 pointN;
//	convertOpenGLCoordinatesToWorldPosition(viewport, projection, modelview,x,y,0,&pointN);
//	ray.first=pointN;
//	ray.second=point;
//}

//-----------------------------------------------------------------------------
// Add connected component for the skeleton construction
//
void SVGshader::ccbForSkeleton (Arrangement_2::Ccb_halfedge_const_circulator circ, SsBuilder& ssb)
{
    Arrangement_2::Ccb_halfedge_const_circulator  curr = circ;
    Arrangement_2::Halfedge_const_handle          he;
    std::vector<Kernel2::Point_2> contourForSsb;
    do
    {
        he = curr;
        if(he->source()->point() == *(he->curve().begin())){
            Curve_2::const_iterator         pit;
            for(pit = he->curve().begin(); pit != he->curve().end(); pit++){
                double x = CGAL::to_double( pit->x());
                double y = CGAL::to_double( pit->y());
                Kernel2::Point_2 pt(x,y);
                contourForSsb.push_back(pt);
            }
        }
        else
        {
            Curve_2::const_reverse_iterator         pit;
            for(pit = he->curve().rbegin(); pit != he->curve().rend(); pit++){
                double x = CGAL::to_double( pit->x());
                double y = CGAL::to_double( pit->y());
                Kernel2::Point_2 pt(x,y);
                contourForSsb.push_back(pt);
            }
        }
        ++curr;
    } while (curr != circ);
    contourForSsb.pop_back();
    try
    {
        ssb.enter_contour(contourForSsb.begin(), contourForSsb.end());
    }
    catch( std::exception const& e )
    {
        std::cout << "EXCEPTION THROWN (" << e.what() << ") during entering contour in ss construction.";
    }

}

//the DC is just for DEBUGGING, it will be deleted later!
void SVGshader::classifyFaces(wxDC & dc, std::vector<SVGshader::FaceType>& faceTypes)
{
    if(myPs->size() == 0 ) return;
    if(arr.number_of_faces() == 0) return;

    std::vector<std::vector< std::pair< gmVector3, gmVector3 > > >::const_iterator iterFaceRaySets;
    std::vector< std::pair<gmVector3, gmVector3> > ::const_iterator iterFaceRays;

    std::vector<std::vector< Kernel2::Point_2> >::const_iterator tempIterVpoints;
    std::vector< Kernel2::Point_2 > ::const_iterator tempIterFpoints;



    std::vector<Implicit*> implicits;
    findImplicits(&implicits);
    if(implicits.size() == 0) return;



    ////To use particle visibility was a BIG HACK, because we cannot simply take any bounding box!
    ////for the visibility of the particles themselves we want small ones...
    ////for the classification we need an object bounding box!
    ////The only thing I see to make this really acceptable would be to use the glMatrix and recover the position...
    ////let's try this for the moment and then we try to solve the global problem later...
    //finding attribute for ParticleBoundingBoxIntersection
    //ParticleVisibility *pVis=0;
    //for (unsigned int s=0;s<myPs->size();++s)
    //{
    //	for (unsigned int p=0;p<(*myPs)[s]->particles.size();++p)
    //	{
    //
    //		pVis= (*myPs)[s]->particles[p]->getAttributeByType<ParticleVisibility>();
    //		if (pVis)
    //			break;
    //	}
    //	if (pVis)
    //		break;
    //}
    //if(!pVis) return;

    //loop over faces
    Arrangement_2::Face_const_iterator    fit;
    for (	fit = arr.faces_begin(),
            tempIterVpoints = interiorPoints.begin(),
            iterFaceRaySets = classifyFaceRaySets.begin()
            ;
            fit != arr.faces_end()
            ;
            ++fit, ++tempIterVpoints, ++iterFaceRaySets)
    {
        if (fit->is_unbounded())
        {
            faceTypes.push_back(FaceType_HOLE);
            continue;
        }

        int isIntersect = 0;
        int isShadow = 0;
        int isSpecularity= 0;

        FaceType ft = FaceType_VISIBLE;
        //if no skeleton exists in this face, damn it as visible
        const std::vector<std::pair<gmVector3,gmVector3> > & faceRays=(*iterFaceRaySets);
        if(faceRays.size() == 0){
            faceTypes.push_back(ft);
            continue;
        }

        for(	iterFaceRays = faceRays.begin()
//                , tempIterFpoints=tempIterVpoints->begin()
                ;
                iterFaceRays != faceRays.end()
                ;
                iterFaceRays++
//                , tempIterFpoints++
                )
        {
            std::pair<gmVector3, Implicit*> result;
            FindFirstIntersectionWithImplicits::findIntersection(iterFaceRays->first, iterFaceRays->second, implicits, &result,100,50,gmEPSILON,10);
            //if no intersection go to the next point
            if(!result.second)	continue;
            //if intersection found, classify it
            FaceType ft = classifyIntersection(result.first);
            isIntersect++;
            switch ( ft )
            {
            case FaceType_SHADOW:
                isShadow++;
                //dc.SetPen(*wxBLACK_PEN);
                //dc.DrawCircle((int) tempIterFpoints->x(), (int) tempIterFpoints->y(),2 );
                break;
            case FaceType_SPECULARITY:
                isSpecularity++;
                //dc.SetPen(*wxGREEN_PEN);
                //dc.DrawCircle((int) tempIterFpoints->x(), (int) tempIterFpoints->y(),2 );
                break;
            default:
                //dc.SetPen(*wxRED_PEN);
                //dc.DrawCircle((int) tempIterFpoints->x(), (int) tempIterFpoints->y(),2 );
                ;
            }
        }
        //checking proportions of visibility events
#ifdef _DEBUG
        std::cout <<"Intersected percentage:"<<(double)isIntersect/(double)faceRays.size()<<'\n';
        std::cout <<"In shadow percentage:"<<(double)isShadow/(double)faceRays.size()<<'\n';
        std::cout <<"In specularity percentage:"<<(double)isSpecularity/(double)faceRays.size()<<'\n';
        std::cout <<"________________________________"<<'\n';
#endif	
        if((double)isIntersect/(double)faceRays.size() < 0.4){
            ft = FaceType_HOLE;
            faceTypes.push_back(ft);
            continue;
        }
        if((double)isShadow/(double)faceRays.size() > 0.6){
            ft = FaceType_SHADOW;
            faceTypes.push_back(ft);
            continue;
        }
        if((double)isSpecularity/(double)faceRays.size() > 0.6){
            ft = FaceType_SPECULARITY;
            faceTypes.push_back(ft);
            continue;
        }

        //no shadow, no spec but still visible, must be lit body
        ft = FaceType_VISIBLE;
        faceTypes.push_back(ft);
    }//end iteration over faces
}

SVGshader::FaceType SVGshader::classifyIntersection(const gmVector3& point)
{
    //test shadow
    if(behaviorShadow){
        if(behaviorShadow->definingManifoldProc(point)<=gmEPSILON) return FaceType_SHADOW;
    }
    //test specularity
    if(behaviorSpecular){
        if (behaviorSpecular->definingManifoldProc(point)<=gmEPSILON) return FaceType_SPECULARITY;
    }
    return FaceType_VISIBLE;
}


void SVGshader::drawSVGnew(wxDC& dc){

    if(!myPs) return;
    if(myPs->size() == 0 ) return;
    if(arr.number_of_faces() == 0) return;
    if (classifyFaceRaySets.size()==0) return;

    //initialize the canvas
    dc.SetBackground(*wxWHITE_BRUSH);

    //set brushes
    wxBrush wxBvisible, wxBshadow, wxBspecularity, wxBhole;
    wxBhole.SetColour(200, 200, 200);
    //set visible brush
    if(shaderSilhouette && shaderSilhouette->IsShaded()){
        std::vector<short> shadeSilhColor;
        shaderSilhouette->getBrushColor(shadeSilhColor);
        wxBvisible.SetColour(shadeSilhColor[0],shadeSilhColor[1],shadeSilhColor[2]);
    }

    //set shadow brush
    if(shaderShadow && shaderShadow->IsShaded()){
        std::vector<short> shadeColor;
        shaderShadow->getShadeColor(shadeColor);
        wxBshadow.SetColour(shadeColor[0],shadeColor[1],shadeColor[2]);
    }

    //set specular brush
    if(shaderSpecular && shaderSpecular->IsShaded()){
        std::vector<short> shadeColor;
        shaderSpecular->getBrushColor(shadeColor);
        wxBspecularity.SetColour(shadeColor[0],shadeColor[1],shadeColor[2]);
    }

    std::vector<SVGshader::FaceType> faceTypes;
    std::vector<SVGshader::FaceType>::const_iterator iterFaceTypes;
    classifyFaces(dc, faceTypes);
    if(faceTypes.size() == 0) return;

    int nT = faceTypes.size();
    int nF = arr.number_of_faces();
    if(nT != nF) return;

    //loop over faces and draw
    Arrangement_2::Face_const_iterator    fit;


    for (fit = arr.faces_begin(), iterFaceTypes = faceTypes.begin(); fit != arr.faces_end(); ++fit, ++iterFaceTypes)
    {
        if (fit->is_unbounded()) continue;
        switch ( *iterFaceTypes )
        {
        case FaceType_HOLE:
            draw_face(fit, dc, wxBhole, false);
            break;
        case FaceType_VISIBLE:
            draw_face(fit, dc, wxBvisible, false);
            break;
        case FaceType_SHADOW:
            draw_face(fit, dc, wxBshadow, false);
            break;
        case FaceType_SPECULARITY:
            draw_face(fit, dc, wxBspecularity, false);
            break;
        default:
            ;
        }
    }


    //draw silhouettes
    wxPen wP;
    if(drawSilhouettesFlag && shaderSilhouette){
        //set pen
        std::vector<short> penColor;
        shaderSilhouette->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderSilhouette->getPenWidth());
        drawSVGContours(dc, wP, myTRANSPARENT_BRUSH, shaderSilhouette, silPts, true);
    }
    //draw shadows
    if(drawShadowsFlag && shaderShadow){
        //set pen
        std::vector<short> penColor;
        shaderShadow->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderShadow->getPenWidth());
        //draw contours
        drawSVGContours(dc, wP, myTRANSPARENT_BRUSH, shaderShadow, shadPts, true);
    }
    //draw suggestive contours
    if(drawSuggestivesFlag && shaderSuggestive){
        //set pen
        std::vector<short> penColor;
        shaderSuggestive->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderSuggestive->getPenWidth());
        //draw contours
        drawSVGContours(dc, wP, myTRANSPARENT_BRUSH, shaderSuggestive, suggPts, true);
    }

    //draw expensive suggestive contours
    if(drawDeCarloSuggestivesFlag && shaderDeCarloSuggestive){
        //set pen
        std::vector<short> penColor;
        shaderDeCarloSuggestive->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderDeCarloSuggestive->getPenWidth());
        //draw contours
        drawSVGContours(dc, wP, myTRANSPARENT_BRUSH, shaderDeCarloSuggestive, deCarloSuggPts, true);
    }

    //draw parabolic contours
    if(drawParabolicsFlag && shaderParabolic){
        //set pen
        std::vector<short> penColor;
        shaderParabolic->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderParabolic->getPenWidth());
        //draw contours
        drawSVGContours(dc, wP, myTRANSPARENT_BRUSH, shaderParabolic, parabPts, true);
    }
    //draw specular highlights
    if(drawSpecularsFlag && shaderSpecular){
        //set pen
        std::vector<short> penColor;
        shaderSpecular->getContourColor(penColor);
        wP.SetColour(penColor[0],penColor[1],penColor[2]);
        wP.SetWidth(shaderSpecular->getPenWidth());
        drawSVGContours(dc, wP, myTRANSPARENT_BRUSH, shaderSpecular, specPts, true);
    }

}

void SVGshader::setSilhouettesFlag(const bool b)
{
    drawSilhouettesFlag=b;
}

void SVGshader::setShadowsFlag(const bool b)
{
    drawShadowsFlag=b;
}
void SVGshader::setSpecularsFlag(const bool b)
{
    drawSpecularsFlag=b;
}
void SVGshader::setSuggestivesFlag(const bool b)
{
    drawSuggestivesFlag=b;
}
void SVGshader::setDeCarloSuggestivesFlag(const bool b)
{
    drawDeCarloSuggestivesFlag=b;
}
void SVGshader::setParabolicsFlag(const bool b)
{
    drawParabolicsFlag=b;
}


#endif
