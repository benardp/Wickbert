#ifndef DISPLAYWINDOW_H
#define DISPLAYWINDOW_H

#include "wx/wx.h"
#include "wx/glcanvas.h"
#include <utility>

#ifdef WB_USE_CGAL
class SVGshader;
#endif


class Modeler;	/* The app class */
class Particles;

/** \class DisplayWindow
 * This class implements an OpenGL viewer intended to display (and interact
 * with) the particle system.
 */
class DisplayWindow : public wxGLCanvas {
public:
	float xrot;
	float yrot;
	float zoom;
	
	Particles *selected_p;
	int selected_particle;

	DisplayWindow(wxWindow *parent, wxWindowID id = wxID_ANY,
        	const wxPoint& pos = wxDefaultPosition,
        	const wxSize& size = wxDefaultSize, long style = 0,
        	const wxString& name = wxT("DisplayWindow"));
	~DisplayWindow();
	void render();
	bool pick(int x, int y);
	void resize();
	void SetupView(int picking, int x, int y);

	void OnPaint(wxPaintEvent& event);
	void OnSize(wxSizeEvent& event);
	void OnEraseBackground(wxEraseEvent& event);
	void OnMouseEvent(wxMouseEvent& event);
	void OnKey(wxKeyEvent& event);
//	void OnKeyDown(wxKeyEvent& event);
	
	virtual void OnDraw(wxDC& dc);
	void OnDrawPlanarMap(wxDC& dc);
	void BuildSkeleton(wxDC& dc);
	void drawSVGnew(wxDC& dc);
	void renderConstrainedDelaunay(wxDC& dc);
	void CaptureGL(const wxString&);

#ifdef WB_USE_CGAL
	SVGshader *svgShader;
#endif

private:
	
	void captureSVG();
	void captureCanvas();
	void initGL();	
	int frame_count;
	bool do_init;

    wxGLContext* glContext;

		
DECLARE_EVENT_TABLE()
};

#endif
