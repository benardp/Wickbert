/**
* Implementation of the Log frame
 * @file LogFrame.cpp
 * @date 08/08/2006
 * @author Matei N. Stroila
 * @remarks
 */

#include "LogFrame.h"
#include "wxModeler.h"

BEGIN_EVENT_TABLE(LogFrame, wxFrame)
	EVT_MENU(wxID_OPEN, LogFrame::OnOpen)
	EVT_MENU(wxID_SAVE, LogFrame::OnSave)
	EVT_MENU(ID_CLEAR, LogFrame::OnClear)
	EVT_MENU(ID_CLOSEWINDOW, LogFrame::OnCloseWindow)
	EVT_CLOSE(LogFrame::OnClose)
END_EVENT_TABLE()

///////////////////////////////////////////////
//LogFrame Implementation (a log window)
///////////////////////////////////////////////

LogFrame::LogFrame(const wxString & title, const wxPoint & pos, const wxSize & size)
: wxFrame(NULL, -1, title, pos, size)
{

	wxMenuBar *menubar = new wxMenuBar;
	
	wxMenu *filemenu = new wxMenu;
	filemenu->Append(wxID_OPEN, "&Open...");
	filemenu->Append(wxID_SAVE, "&Save...");
	filemenu->Append(ID_CLEAR, "&Clear");
	filemenu->Append(ID_CLOSEWINDOW, "&Close");
	
	menubar->Append(filemenu,"&File");
	
	SetMenuBar(menubar);
	
	logTxtControl = new wxTextCtrl(this, -1, "", wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE);
	redirect=new wxStreamToTextRedirector(logTxtControl);
}

void LogFrame::OnOpen(wxCommandEvent &event)
{
	
	wxFileDialog *filedialog = new wxFileDialog(this, "Open a File", "", "",
												"Txt (*.txt)|*.txt|"
												"Other files (*.*)|*.*", wxFD_OPEN);
	
	if (filedialog->ShowModal() == wxID_CANCEL)
		return;
	
	wxString filename = filedialog->GetPath();
	logTxtControl->LoadFile(filename);
	
}

void LogFrame::OnSave(wxCommandEvent &event)
{
	
	wxFileDialog *filedialog = new wxFileDialog(this, "Save to a File", "", "",
												"Text (*.txt)|*.txt|"
												"Other files (*.*)|*.*", wxFD_SAVE);
	
	if (filedialog->ShowModal() == wxID_CANCEL)
		return;
	
	wxString filename = filedialog->GetPath();
	
	std::ofstream out(filename.ToStdString());
	if (out) {
		out << logTxtControl->GetValue() << std::endl;
		out.close();	
	}
}

void LogFrame::OnClear(wxCommandEvent &event)
{
	logTxtControl->Clear();
}

void LogFrame::OnCloseWindow(wxCommandEvent& event)
{
	wxMessageDialog* dialog = new wxMessageDialog(this,
												  "Save log data?", "wxModeler", wxYES_NO|wxCANCEL);
	
	int ans = dialog->ShowModal();
	dialog->Destroy();
	
	switch (ans){
        case wxID_YES:      // Save, then destroy, quitting app
			OnSave(event);
			this->Destroy();
			wxGetApp().logWindow = NULL; 
			break;
        case wxID_NO:       // Don't save; just destroy, quitting app
			this->Destroy();
			wxGetApp().logWindow = NULL; 
			break;
        case wxID_CANCEL:   // Do nothing - so don't quit app.
        default:
			break;
	}
	delete this->redirect;
}


void LogFrame::OnClose(wxCloseEvent &event)
{
	this->Destroy();
	delete this->redirect;
	wxGetApp().logWindow = NULL; 
}

