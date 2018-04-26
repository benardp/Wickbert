/**
* Declaration of the Log frame
 * @file LogFrame.h
 * @date 08/08/2006
 * @author Matei N. Stroila
 * @remarks
 */

#ifndef LOGFRAME_H
#define LOGFRAME_H

#include "wx/wx.h"

/** \class LogFrame
* This is a top-level frame for logging stuff in a text ctrl.
*/

class LogFrame : public wxFrame
{
protected:
	wxStreamToTextRedirector * redirect;

public:
	
	LogFrame(const wxString& title, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize);
	void OnClose(wxCloseEvent &event);
	void OnOpen(wxCommandEvent& event);
	void OnSave(wxCommandEvent& event);
	void OnClear(wxCommandEvent& event);
	void OnCloseWindow(wxCommandEvent& event);
	
	wxTextCtrl *logTxtControl; //text control that goes into the window log --ms 
	
	DECLARE_EVENT_TABLE()
};

#endif
