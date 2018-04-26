/**
 * Declaration of an SVG shader 
 * @file SVGshader.h
 * @date 04/08/2006
 * @author Matei N. Stroila
 * @remarks
 */

#ifdef WB_USE_CGAL
#ifndef SVGSHADER_H
#define SVGSHADER_H

//#define CGAL_MP_FLOAT_ALLOW_INEXACT
//#define CGAL_NO_PRECONDITIONS

//CGAL headers
#include <CGAL/basic.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>

//Skeleton headers
#include<CGAL/Polygon_2.h>
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Straight_skeleton_builder_2.h>

//Constrained Delaunay triangulation headers
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/number_utils_classes.h>
#include <CGAL/circulator.h>

#include "wx/wx.h"
#include "wx/dcsvg.h" 
#include "Particles/Attributes/Contours.h"

#include <boost/shared_ptr.hpp>



// Kernel for Arrangments
typedef CGAL::Gmpq               Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef Kernel::FT										NT;
typedef CGAL::Arr_segment_traits_2<Kernel>              Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>   Polyline_traits_2;
typedef Polyline_traits_2::Curve_2                      Polyline_2;

typedef CGAL::Arr_consolidated_curve_data_traits_2<Polyline_traits_2, int> Traits_2;

typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Curve_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2>     PointLocation;


//skeleton
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel2;
typedef CGAL::Polygon_2<Kernel2>    Contour;
typedef boost::shared_ptr<Contour> ContourPtr;
typedef std::vector<ContourPtr>    ContourSequence ;
typedef CGAL::Straight_skeleton_2<Kernel2> Ss;
typedef CGAL::Straight_skeleton_builder_traits_2<Kernel2>      SsBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;

//Constrained Delaunay triangulation
typedef CGAL::Triangulation_vertex_base_2<Kernel>                     Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel>           Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag> CDT;

class ShaderSilhouetteContour;
class SilhouetteAdhesion;
class ShaderShadowContour;
class ShadowAdhesion;
class ShaderSpecularContour;
class SpecularAdhesion;
class ShaderSuggestiveContour;
class SuggestiveContourAdhesion;
class ShaderParabolicsContour;
class ParabolicsAdhesion;
class DeCarloSuggestiveContourAdhesion;
class ShaderDeCarloSuggestiveContour;
class ShaderContour;
class Implicit;

class SVGshader 
{

	typedef std::vector<Point_2> VectorOfPoint_2;
	typedef std::pair< int, int> ContourEnd; //loop index - particle index
	typedef std::pair<ContourEnd, ContourEnd> EndsOfContour; 

	enum FaceType {FaceType_HOLE, FaceType_VISIBLE, FaceType_SHADOW, FaceType_SPECULARITY}; 

public:
	
	SVGshader();
    ~SVGshader(){}

	void setParticleSystems(ParticleSystems*);
	void drawSVG(wxDC& dc);	
	void builAndDrawPlanarMap(wxDC& dc);
	void buildSkeleton(wxDC& dc);
	void drawSVGnew(wxDC& dc);
	void buildAndDrawConstrainedDelaunay(wxDC& dc);

	void setSilhouettesFlag(const bool);
	void setShadowsFlag(const bool);
	void setSpecularsFlag(const bool);
	void setSuggestivesFlag(const bool);
	void setDeCarloSuggestivesFlag(const bool);
	void setParabolicsFlag(const bool);

private:
	//later these two functions should be IN wxDC.

	ParticleSystems *myPs;

	std::vector<Point_2> silPts;
	std::vector<Point_2> shadPts;
	std::vector<Point_2> specPts;
	std::vector<Point_2> suggPts;
	std::vector<Point_2> deCarloSuggPts;
	std::vector<Point_2> parabPts;
	
	//in the next step we should not really store these, but just a list of Shader/Adhesion couples
	const ShaderSilhouetteContour* shaderSilhouette;
	const SilhouetteAdhesion* behaviorSilhouette;
	
	const ShaderShadowContour* shaderShadow;
	const ShadowAdhesion * behaviorShadow;

	const ShaderSpecularContour* shaderSpecular;
	const SpecularAdhesion * behaviorSpecular;

	const ShaderSuggestiveContour* shaderSuggestive;
	const SuggestiveContourAdhesion * behaviorSuggestive;

	const ShaderParabolicsContour* shaderParabolic;
	const ParabolicsAdhesion * behaviorParabolic;

	const ShaderDeCarloSuggestiveContour* shaderDeCarloSuggestive;
	const DeCarloSuggestiveContourAdhesion * behaviorSuggestiveExp;

	bool projectLoopsOnNormalizedCoords(const ShaderContour* shader, std::vector<Point_2>& pts);
	void findContoursAndBehaviors();
	void projectContoursOnNormalizedCoords();

	

	void drawCycles(wxDC& dc, 
							const ShaderContour* shaderSil, 
							const ShaderContour* shaderContour,
							const std::vector<Point_2>& contPts);
	
	void computeVisibleContours(const ShaderContour* shaderSil, 
								const ShaderContour* shaderContour,
								const std::vector<Point_2> contour2DPts,
								std::vector<VectorOfPoint_2>& visContours,
								std::vector<EndsOfContour>& contourEnds);

	bool findClosestParticle(const Point_2& sourcePt, 
							const ShaderContour* shaderSil, 
							Point_2& targetPt,
							int& p_index, //the closest particle index
							int& l_index); //the corresponding loop 
	
	/** 
	* Find all implicits in the scene	
	* @param implicits a vector of implicit objects. Not cleared by this function.
	*/	
	void findImplicits(std::vector<Implicit*>* implicits);
		
	Arrangement_2   arr;
	//world coordinates of rays that we want to test to classify the faces
	std::vector<std::vector< std::pair<gmVector3, gmVector3> > > classifyFaceRaySets;
	//some debug display.
	std::vector<std::vector<Kernel2::Point_2> > interiorPoints;
	
	void drawSVGContours(wxDC& dc, 
						 wxPen& pen, 
						 wxBrush& brush,  
						 const ShaderContour*, 
						 const std::vector<Point_2>& contPts,
						 const bool drawContoursOnly=false);
	
	unsigned int addFaceVerticesInDCCoordinates(const wxDC& dc, Arrangement_2::Ccb_halfedge_const_circulator circ, std::vector<wxPoint>& points);
	void draw_face (Arrangement_2::Face_const_handle f, wxDC& dc, 
					wxBrush& wB, 
					const bool randomColour);
	void draw_arrangement (wxDC& dc);
	
	void constructSkeleton (Arrangement_2::Face_const_handle f, wxDC& dc);
	void ccbForSkeleton (Arrangement_2::Ccb_halfedge_const_circulator circ, SsBuilder& ssb);

	void classifyFaces(wxDC& dc, std::vector<FaceType>& faceTypes);
	FaceType classifyIntersection(const gmVector3& point);

	//Constrained Delaunay triangulation
	CDT cdt;

	//void computeFaceRay(std::pair<gmVector3, gmVector3>& ray, wxDC& dc, const double vx, const double vy, 
	//				 double *projection,   double *modelview,   int *viewport );


	bool drawSilhouettesFlag;
	bool drawShadowsFlag;
	bool drawSpecularsFlag;
	bool drawSuggestivesFlag;
	bool drawDeCarloSuggestivesFlag;
	bool drawParabolicsFlag;
};

#endif
#endif

