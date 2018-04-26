/////////////////////////////////////////////////////////////////////////////
// Viewer is based on...
// Name:        pngdemo.cpp
// Purpose:     Demos PNG reading
// Author:      Julian Smart
// Modified by:
// Created:     04/01/98
// RCS-ID:      $Id: pngdemo.cpp,v 1.17.2.1 2002/12/15 17:25:21 MBN Exp $
// Copyright:   (c) Julian Smart and Markus Holzem
// Licence:     wxWindows license
/////////////////////////////////////////////////////////////////////////////

#ifdef __GNUG__
#pragma implementation "viewer.h"
#endif

// For compilers that support precompilation, includes "wx.h".
//#include "wx/wxprec.h"
#include "wx/wx.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#include "wx/image.h"

#include "viewer.h"

MyFrame   *frame = (MyFrame *) NULL;
wxBitmap  *g_TestBitmap = (wxBitmap *) NULL;

IMPLEMENT_APP(MyApp)

MyApp::MyApp()
{
}

bool MyApp::OnInit(void)
{
  wxImage::AddHandler(new wxPNGHandler);

  // Create the main frame window
  frame = new MyFrame((wxFrame *) NULL, _T("viewer"), wxPoint(0, 0), wxSize(300, 300));

  // Give it a status line
  frame->CreateStatusBar(2);

  // Make a menubar
  wxMenu *file_menu = new wxMenu;
  wxMenu *help_menu = new wxMenu;

  file_menu->Append(VIEWER_LOAD_FILE, _T("&Load file"),                _T("Load file"));
  file_menu->Append(VIEWER_SAVE_FILE, _T("&Save file"),                _T("Save file"));
  file_menu->Append(VIEWER_QUIT, _T("E&xit"),                _T("Quit program"));
  help_menu->Append(VIEWER_ABOUT, _T("&About"),              _T("About viewer"));

  wxMenuBar *menu_bar = new wxMenuBar;

  menu_bar->Append(file_menu, _T("&File"));
  menu_bar->Append(help_menu, _T("&Help"));

  // Associate the menu bar with the frame
  frame->SetMenuBar(menu_bar);

  MyCanvas *canvas = new MyCanvas(frame, wxPoint(0, 0), wxSize(100, 100));

  // Give it scrollbars: the virtual canvas is 20 * 50 = 1000 pixels in each direction
//  canvas->SetScrollbars(20, 20, 50, 50, 4, 4);
  frame->canvas = canvas;

  frame->Show(TRUE);

  frame->SetStatusText(_T("Convert BMP to XPM"));

  return TRUE;
}

BEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_MENU(VIEWER_QUIT,      MyFrame::OnQuit)
    EVT_MENU(VIEWER_ABOUT,     MyFrame::OnAbout)
    EVT_MENU(VIEWER_LOAD_FILE, MyFrame::OnLoadFile)
    EVT_MENU(VIEWER_SAVE_FILE, MyFrame::OnSaveFile)
END_EVENT_TABLE()

// Define my frame constructor
MyFrame::MyFrame(wxFrame *frame, const wxString& title, const wxPoint& pos, const wxSize& size):
  wxFrame(frame, -1, title, pos, size)
{
  canvas = (MyCanvas *) NULL;
}

// frame destructor
MyFrame::~MyFrame()
{
    if (g_TestBitmap)
    {
        delete g_TestBitmap;
        g_TestBitmap = (wxBitmap *) NULL;
    }
}

void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
    Close(TRUE);
}

void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
    (void)wxMessageBox(_T("Adapted from\nPNG demo\nJulian Smart (c) 1998"),
            _T("About viewer"), wxOK);
}

void MyFrame::OnSaveFile(wxCommandEvent& WXUNUSED(event))
{
#if 0
  wxString f = wxFileSelector( wxT("Save Image"), (const wxChar *)NULL,
                               (const wxChar *)NULL,
                               wxT("xpm"), wxT("XPM files (*.xpm)|*.xpm") );
  if (f == _T(""))  return;
#else
	wxFileDialog *fd = new wxFileDialog(this, _T("Save as XPM file"), "", "", "*.xpm", wxSAVE);
	if (fd->ShowModal() != wxID_OK) return;
	wxString f = fd->GetPath();
#endif


#if 0

  wxBitmap *backstore = new wxBitmap( 150, 150 );
  
  wxMemoryDC memDC;
  memDC.SelectObject( *backstore );
  memDC.Clear();
  memDC.SetBrush( *wxBLACK_BRUSH );
  memDC.SetPen( *wxWHITE_PEN );
  memDC.DrawRectangle( 0, 0, 150, 150 );
  memDC.SetPen( *wxBLACK_PEN );
  memDC.DrawLine( 0, 0, 0, 10 );
  memDC.SetTextForeground( *wxWHITE );
  memDC.DrawText( _T("This is a memory dc."), 10, 10 );
  
  memDC.SelectObject( wxNullBitmap );
  
  backstore->SaveFile( f, wxBITMAP_TYPE_XPM, (wxPalette*)NULL );
  
  delete backstore;
#else
  g_TestBitmap->SaveFile(f, wxBITMAP_TYPE_XPM, (wxPalette*)NULL );
#endif
}

void MyFrame::OnLoadFile(wxCommandEvent& WXUNUSED(event))
{
    // Show file selector.
    wxString f = wxFileSelector(wxT("Open Image"), (const wxChar *) NULL,
                                    (const wxChar *) NULL, wxT("bmp"),
                                    wxT("BMP files (*.bmp)|*.bmp"));

    if (f == _T(""))
        return;

    if ( g_TestBitmap )
        delete g_TestBitmap;

	g_TestBitmap = new wxBitmap(f, wxBITMAP_TYPE_BMP);
    if (!g_TestBitmap->Ok())
    {
        delete g_TestBitmap;
        g_TestBitmap = (wxBitmap *) NULL;
    }

    canvas->Refresh();
}

BEGIN_EVENT_TABLE(MyCanvas, wxScrolledWindow)
    EVT_PAINT(MyCanvas::OnPaint)
END_EVENT_TABLE()

// Define a constructor for my canvas
MyCanvas::MyCanvas(wxWindow *parent, const wxPoint& pos, const wxSize& size):
 wxScrolledWindow(parent, -1, pos, size)
{
}

MyCanvas::~MyCanvas(void)
{
}

// Define the repainting behaviour
void MyCanvas::OnPaint(wxPaintEvent& WXUNUSED(event))
{
  wxPaintDC dc(this);
  dc.SetPen(* wxRED_PEN);

  int i;
  for ( i = 0; i < 500; i += 10)
  {
    dc.DrawLine(0, i, 800, i);
  }
  if ( g_TestBitmap && g_TestBitmap->Ok() )
  {
    wxMemoryDC memDC;
    if ( g_TestBitmap->GetPalette() )
    {
        memDC.SetPalette(* g_TestBitmap->GetPalette());
        dc.SetPalette(* g_TestBitmap->GetPalette());
    }
    memDC.SelectObject(* g_TestBitmap);

    // Normal, non-transparent blitting
    dc.Blit(20, 20, g_TestBitmap->GetWidth(), g_TestBitmap->GetHeight(), & memDC, 0, 0, wxCOPY, FALSE);

    memDC.SelectObject(wxNullBitmap);
  }

  if ( g_TestBitmap && g_TestBitmap->Ok() )
  {
    wxMemoryDC memDC;
    memDC.SelectObject(* g_TestBitmap);

    // Transparent blitting if there's a mask in the bitmap
    dc.Blit(20 + g_TestBitmap->GetWidth() + 20, 20, g_TestBitmap->GetWidth(), g_TestBitmap->GetHeight(), & memDC,
      0, 0, wxCOPY, TRUE);

    memDC.SelectObject(wxNullBitmap);
  }
}


