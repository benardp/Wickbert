/**
 * Implementation of the Clipart frame
 * @file ClipartFrame.cpp
 * @date 08/08/2006
 * @author Matei N. Stroila
 * @remarks
 */

#ifdef WB_USE_CGAL

#include "ClipartFrame.h"
#include "wxModeler.h"
#include "DisplayWindow.h"
#include "SVGshader.h"

#include "bitmaps/exefile.xpm"
#include "bitmaps/filesave.xpm"
#include "wx/dcbuffer.h"

BEGIN_EVENT_TABLE(ClipArtFrame, wxFrame)
	EVT_MENU(wxID_SAVE, ClipArtFrame::OnSave)
	EVT_MENU(ID_CLOSE_CLIPART_WINDOW, ClipArtFrame::OnCloseWindow)
	EVT_MENU(ID_CLIPART_DRAW_SILHOUETTES, ClipArtFrame::OnDrawSilhouettes)
	EVT_MENU(ID_CLIPART_DRAW_SHADOWS, ClipArtFrame::OnDrawShadows)
	EVT_MENU(ID_CLIPART_DRAW_SPECULARS, ClipArtFrame::OnDrawSpeculars)
	EVT_MENU(ID_CLIPART_DRAW_SUGGESTIVES, ClipArtFrame::OnDrawSuggestives)
	EVT_MENU(ID_CLIPART_DRAW_DECARLO_SUGGESTIVES, ClipArtFrame::OnDrawDeCarloSuggestives)
	EVT_MENU(ID_CLIPART_DRAW_PARABOLICS, ClipArtFrame::OnDrawParabolics)

	EVT_CLOSE(ClipArtFrame::OnClose)
	EVT_TOOL(ID_RENDER, ClipArtFrame::OnRender)
	EVT_TOOL(ID_RENDER_PLANARMAP, ClipArtFrame::OnRenderPlanarMap)
	EVT_TOOL(ID_BUILD_SKELETON, ClipArtFrame::OnBuildSkeleton)
	EVT_TOOL(ID_RENDER_SVG, ClipArtFrame::OnRenderSVG)
	EVT_TOOL(ID_RENDER_CDT, ClipArtFrame::OnRenderConstrainedDelaunay)
	EVT_ERASE_BACKGROUND(ClipArtFrame::OnErase)
END_EVENT_TABLE()


ClipArtFrame::ClipArtFrame(const wxString &title, const wxPoint &pos, const wxSize &size)
: wxFrame(NULL, -1, title, pos, size)
{
	wxMenuBar *menubar = new wxMenuBar;
	
	wxMenu *filemenu = new wxMenu;
	filemenu->Append(wxID_SAVE, "&Save...");
	filemenu->Append(ID_CLOSE_CLIPART_WINDOW, "&Close");

	wxMenu *drawMenu = new wxMenu;
	drawMenu->AppendCheckItem(ID_CLIPART_DRAW_SILHOUETTES, "Draw silhouettes");
	drawMenu->AppendCheckItem(ID_CLIPART_DRAW_SHADOWS, "Draw shadows");
	drawMenu->AppendCheckItem(ID_CLIPART_DRAW_SPECULARS, "Draw speculars");
	drawMenu->AppendCheckItem(ID_CLIPART_DRAW_SUGGESTIVES, "Draw suggestives");
	drawMenu->AppendCheckItem(ID_CLIPART_DRAW_DECARLO_SUGGESTIVES, "Draw DeCarlo suggestives");
	drawMenu->AppendCheckItem(ID_CLIPART_DRAW_PARABOLICS, "Draw parabolics");

	
	menubar->Append(filemenu,"&File");
	menubar->Append(drawMenu,"&Drawing");
	
	SetMenuBar(menubar);
	
	// Create the toolbar
	this->CreateToolBar(wxNO_BORDER|wxHORIZONTAL|wxTB_FLAT, ID_TOOLBAR);
	
	this->GetToolBar()->SetMargins( 2, 2 );
	
	this->InitToolbar(this->GetToolBar());
	
	
}

bool ClipArtFrame::InitToolbar(wxToolBar* toolBar)
{
	// Set up toolbar
	wxBitmap* toolBarBitmaps[6];
	
	toolBarBitmaps[0] = new wxBitmap( exefile_xpm );
	toolBarBitmaps[1] = new wxBitmap( exefile_xpm );
	toolBarBitmaps[2] = new wxBitmap( exefile_xpm );
	toolBarBitmaps[3] = new wxBitmap( exefile_xpm );
	toolBarBitmaps[4] = new wxBitmap( exefile_xpm );
	toolBarBitmaps[5] = new wxBitmap( filesave_xpm );
	
    toolBar->AddTool(ID_RENDER, "Render", *(toolBarBitmaps[0]), "Render Clipart");
    toolBar->AddTool(ID_RENDER_PLANARMAP, "Planar Map", *(toolBarBitmaps[1]), "Build and Render Planar Map");
    toolBar->AddTool(ID_RENDER_CDT, "CDT", *(toolBarBitmaps[3]), "Build and Render Constrained Delaunay");
    toolBar->AddTool(ID_RENDER_SVG, "SVG", *(toolBarBitmaps[4]), "Render Correct SVG");
    toolBar->AddTool(wxID_SAVE, "Save", *(toolBarBitmaps[5]), "Save file");
	
	
	toolBar->Realize();
	
	// Can delete the bitmaps since they're reference counted
	int i;
	for (i = 0; i < 6; i++)
		delete toolBarBitmaps[i];
	
	return true;
}

void ClipArtFrame::OnRender(wxCommandEvent& event)
{
	//update the clipart window
	wxClientDC dcClipArt(this);
	wxBufferedDC buffer(&dcClipArt); 
	(wxGetApp().mainframe)->displaywin->OnDraw(buffer);
}

void ClipArtFrame::OnRenderSVG(wxCommandEvent& event){
	//update the clipart window
	wxClientDC dcClipArt(this);
	wxBufferedDC buffer(&dcClipArt); 
	(wxGetApp().mainframe)->displaywin->drawSVGnew(buffer);
}

void ClipArtFrame::OnRenderPlanarMap(wxCommandEvent& event)
{
	//update the clipart window
	wxClientDC dcClipArt(this);
	wxBufferedDC buffer(&dcClipArt); 
	(wxGetApp().mainframe)->displaywin->OnDrawPlanarMap(buffer);
}

void  ClipArtFrame::OnRenderConstrainedDelaunay(wxCommandEvent& event)
{
	//update the clipart window
	wxClientDC dcClipArt(this);
	wxBufferedDC buffer(&dcClipArt); 
	(wxGetApp().mainframe)->displaywin->renderConstrainedDelaunay(buffer);
}

void ClipArtFrame::OnBuildSkeleton(wxCommandEvent& event)
{
	//update the clipart window
	wxClientDC dcClipArt(this);
	wxBufferedDC buffer(&dcClipArt); 
	(wxGetApp().mainframe)->displaywin->BuildSkeleton(buffer);
}

void ClipArtFrame::OnErase(wxEraseEvent& event)
{
	//update the clipart window
	//wxClientDC dcClipArt(this);
	//wxBufferedDC buffer(&dcClipArt); 
	//(wxGetApp().mainframe)->displaywin->OnDraw(buffer);
}

void ClipArtFrame::OnSave(wxCommandEvent &event)
{
	
	wxFileDialog *filedialog = new wxFileDialog(this, "Save to a File", "", "", "SVG (*.svg)|*.svg|""Other files (*.*)|*.*", wxFD_SAVE);
	
	if (filedialog->ShowModal() == wxID_CANCEL)
		return;
	wxString filename = filedialog->GetPath();
	
	int w, h;
	(wxGetApp().mainframe)->displaywin->GetClientSize(&w, &h);
	//We virtually increment the resolution
	//the final output will have the approximate size as on the screen (a monitor has about 72 dpi)
	wxSVGFileDC svgDC (filename, 1000*w, 1000*h,1000*72);
	(wxGetApp().mainframe)->displaywin->drawSVGnew(svgDC);
}

void ClipArtFrame::OnCloseWindow(wxCommandEvent& event)
{
	wxMessageDialog* dialog = new wxMessageDialog(this, "Save SVG frame?", "wxModeler", wxYES_NO|wxCANCEL);
	
	int ans = dialog->ShowModal();
	dialog->Destroy();
	
	switch (ans){
        case wxID_YES:      // Save, then destroy, quitting app
			OnSave(event);
			this->Destroy();
			wxGetApp().clipArtFrame = NULL; 
			break;
        case wxID_NO:       // Don't save; just destroy, quitting app
			this->Destroy();
			wxGetApp().clipArtFrame = NULL; 
			break;
        case wxID_CANCEL:   // Do nothing - so don't quit app.
        default:
			break;
	}
}


void ClipArtFrame::OnClose(wxCloseEvent &event)
{
	this->Destroy();
	wxGetApp().clipArtFrame = NULL; 
}

void ClipArtFrame::OnDrawSilhouettes(wxCommandEvent& event)
{
	if((wxGetApp().mainframe)->displaywin->svgShader)
		(wxGetApp().mainframe)->displaywin->svgShader->setSilhouettesFlag(event.IsChecked());
}

void ClipArtFrame::OnDrawShadows(wxCommandEvent& event)
{
	if((wxGetApp().mainframe)->displaywin->svgShader)
		(wxGetApp().mainframe)->displaywin->svgShader->setShadowsFlag(event.IsChecked());
}

void ClipArtFrame::OnDrawSpeculars(wxCommandEvent& event)
{
	if((wxGetApp().mainframe)->displaywin->svgShader)
		(wxGetApp().mainframe)->displaywin->svgShader->setSpecularsFlag(event.IsChecked());
}

void ClipArtFrame::OnDrawSuggestives(wxCommandEvent& event)
{
	if((wxGetApp().mainframe)->displaywin->svgShader)
		(wxGetApp().mainframe)->displaywin->svgShader->setSuggestivesFlag(event.IsChecked());
}

void ClipArtFrame::OnDrawDeCarloSuggestives(wxCommandEvent& event)
{
	if((wxGetApp().mainframe)->displaywin->svgShader)
		(wxGetApp().mainframe)->displaywin->svgShader->setDeCarloSuggestivesFlag(event.IsChecked());
}

void ClipArtFrame::OnDrawParabolics(wxCommandEvent& event)
{
	if((wxGetApp().mainframe)->displaywin->svgShader)
		(wxGetApp().mainframe)->displaywin->svgShader->setParabolicsFlag(event.IsChecked());
}

#endif
