#ifndef WXMODELER_H
#define WXMODELER_H

#include "wx/wx.h"
#include "wx/glcanvas.h"
#include "wx/splitter.h"
#include "wx/tooltip.h"
#include "wx/dcsvg.h"
#include "wx/toolbar.h"

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#ifdef check
#undef check
#endif

#include "Implicit/Implicit.h"		// for the ImplicitRegistry
#include "Particles/ParticleSystem.h"	// For the ParticleStuffRegistry
#include "Surface/Surface.h"

#include <vector>

enum {
	ID_PARTICLETREE = 1,
	ID_PARTICLETREEPARAM,
	ID_PARTICLETREEPARAMBUTTON,
	ID_PARTICLEPARAMCOMBOBOX,
	ID_PARTICLETREEPARAMCHECKBOX,
	ID_PARTICLESREPORT,
	ID_ADDPARTICLE,
	ID_KILLPARTICLES,
	ID_PLAYPARTICLES,
	ID_SHOWPARTICLES,
	ID_SURFACETREE,
	ID_SURFACETREEPARAMBUTTON,
	ID_SURFACETREEPARAMCOMBOBOX,
	ID_IMPLICITPARAM,
	ID_UPDEGREE,
	ID_DOWNDEGREE,
	ID_SURFACEPARAM,
	ID_NEWSURFACE,
	ID_NEWPARTICLESYSTEM,
	ID_NEWPARTICLES,
	ID_COPYPARTICLES,
	ID_NEWPARTICLESTUFF,
	ID_IMPINTCHANGED,
	ID_PAUSEALL,
	ID_PLAYALL,
	ID_REFRESHPARWIN,
	ID_PARPROFILE,
	ID_OPENGLEXTS,
	ID_ZOOMIN,
	ID_ZOOMOUT,
	ID_SAVE_SVG,
	ID_SAVE_GL,
	ID_START_SAVE_SVG_ANIMATION,
	ID_STOP_SAVE_SVG_ANIMATION,
	ID_START_SAVE_GL_ANIMATION,
	ID_STOP_SAVE_GL_ANIMATION,
	ID_CLEAR, //clear log window
	ID_CLOSEWINDOW, //close log window
	ID_CLOSE_CLIPART_WINDOW, //close log window
	ID_QUIT,
	ID_ABOUT,
	ID_TOOLBAR,
	ID_COMBO,
	ID_NEWIMP,
	ID_SHOWLOG,
	ID_SHOWCLIPART,	
	ID_RENDER,
	ID_RENDER_PLANARMAP,
	ID_BUILD_SKELETON,
	ID_RENDER_SVG,
	ID_RENDER_CDT,
	ID_CLIPART_DRAW_SILHOUETTES,
	ID_CLIPART_DRAW_SHADOWS,
	ID_CLIPART_DRAW_SPECULARS,
	ID_CLIPART_DRAW_SUGGESTIVES,
	ID_CLIPART_DRAW_DECARLO_SUGGESTIVES,
	ID_CLIPART_DRAW_PARABOLICS,
	/* These are just placeholders for the new surface menu because wxWindows insists that
	   menu items each have a unique identifier by which they are accessed. */
	ID_NEWSURFACEMENU = 100,
	ID_NEWSURFACEMENU_LAST = 199,
	/* next item should start no less than 200 */

};

class DisplayWindow;
class MainFrame;
class LogFrame;
class ClipArtFrame;
class ParticleWidget;
class SurfaceTree;

/** \class Modeler
 * This is the main application class. It contains the surface and particle
 * system factories, and vectors of surfaces and particle systems.
 */
class Modeler : public wxApp
{
public:
	/** Used to create new implicits
	 */
	ImplicitRegistry *impReg;
	ParticleStuffRegistry *psReg;

	Surfaces surfaces;
	ParticleSystems particlesystems;

	virtual bool OnInit();
	int OnExit();

	MainFrame *mainframe;
	LogFrame *logWindow; //window log --ms
	bool captureSVGflag, captureCanvasFlag; //flag the SVG or GL animation capture --ms

#ifdef WB_USE_CGAL
	ClipArtFrame *clipArtFrame; //clipart window --ms
#endif
	
DECLARE_EVENT_TABLE()
};


DECLARE_APP(Modeler)	/* Declares a Modeler &wxGetApp() function */

/** \class MainFrame
 * This is the main top-level frame of the modeler.
 */
class MainFrame : public wxFrame
{
public:
	DisplayWindow *displaywin;
	ParticleWidget *parwin;
	SurfaceTree *impwin;

	std::vector<std::string *> idtosurfname;

	//for std::cout redirection to log window:
	std::streambuf *cerrbuforig;
	std::streambuf *coutbuforig;

	MainFrame(const wxString& title, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize);

    void OnIdle(wxIdleEvent& event);
	void OnOpen(wxCommandEvent& event);
	void OnSave(wxCommandEvent& event);
	void OnQuit(wxCommandEvent& event);
	void OnAbout(wxCommandEvent& event);
	void OnNewSurface(wxCommandEvent& event);
	void OnNewSurfaceMenu(wxCommandEvent &event);
	void OnPauseAll(wxCommandEvent& event);
	void OnPlayAll(wxCommandEvent& event);
	void OnRefreshParWin(wxCommandEvent& event);
	void OnParProfile(wxCommandEvent& event);
	void OnOpenGLExts(wxCommandEvent& event);
	void OnSaveSVG(wxCommandEvent& event);
	void OnSaveGL(wxCommandEvent& event);
	void OnStartSaveSVGAnimation(wxCommandEvent& event);
	void OnStopSaveSVGAnimation(wxCommandEvent& event);
	void OnStartSaveGLAnimation(wxCommandEvent& event);
	void OnStopSaveGLAnimation(wxCommandEvent& event);
	void OnZoomIn(wxCommandEvent& event);
	void OnZoomOut(wxCommandEvent& event);
	void LoadWB(wxString filename);
	void LoadPar(wxString filename);
	
	void OnShowLog(wxCommandEvent& event); //callback that shows the window log --ms
	void OnShowClipArt(wxCommandEvent& event); //callback that shows the Clip Art window --ms
	
	wxMenu *NewSurfaceMenu();

DECLARE_EVENT_TABLE()
};

#endif

