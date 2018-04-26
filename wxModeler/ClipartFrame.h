/**
 * Declaration of the Clipart frame
 * @file ClipartFrame.h
 * @date 08/08/2006
 * @author Matei N. Stroila
 * @remarks
 */
#ifndef CLIPARTFRAME_H
#define CLIPARTFRAME_H

#include "wx/wx.h"

/** \class ClipArtFrame
 * This is a top-level frame for rendering ClipArt.
 */

class ClipArtFrame : public wxFrame
{
public:
	
	ClipArtFrame(const wxString& title, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize);
	void OnClose(wxCloseEvent &event);
	void OnSave(wxCommandEvent& event);
	void OnCloseWindow(wxCommandEvent& event);
	bool InitToolbar(wxToolBar* toolBar);
	void OnRender(wxCommandEvent& event);
	void OnRenderPlanarMap(wxCommandEvent& event);
	void OnRenderConstrainedDelaunay(wxCommandEvent& event);
	void OnBuildSkeleton(wxCommandEvent& event);
	void OnRenderSVG(wxCommandEvent& event);
	void OnErase(wxEraseEvent& event);
	void OnDrawSilhouettes(wxCommandEvent& event);
	void OnDrawShadows(wxCommandEvent& event);
	void OnDrawSpeculars(wxCommandEvent& event);
	void OnDrawParabolics(wxCommandEvent& event);
	void OnDrawSuggestives(wxCommandEvent& event);
	void OnDrawDeCarloSuggestives(wxCommandEvent& event);

	DECLARE_EVENT_TABLE()
};


#endif