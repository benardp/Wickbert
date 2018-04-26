#ifndef SURFACETREE_H
#define SURFACETREE_H

#include "wx/wx.h"
#include "wx/listctrl.h"
#include "Surface/Implicit/Implicit.h"
#include "Surface/Implicit/Algebraic/Algebraic.h"

class ImpParamTextCtrl : public wxTextCtrl
{
public:
	Implicit *imp;
	int i;

	ImpParamTextCtrl(wxWindow *parent, wxString &value, Implicit *imp, int i)
	: wxTextCtrl(parent,ID_IMPLICITPARAM,value,wxDefaultPosition,wxSize(200,20),wxTE_PROCESS_ENTER)
	{ this->imp = imp; this->i = i; }

	void OnValue(wxCommandEvent &event);
	void OnUpdateUI(wxUpdateUIEvent &event);

DECLARE_EVENT_TABLE()
};

class SurfParamTextCtrl : public wxTextCtrl
{
public:
	Surface *s;
	int i;

	SurfParamTextCtrl(wxWindow *parent, const wxString &value, Surface *s, int i)
	: wxTextCtrl(parent,ID_SURFACEPARAM,value,wxDefaultPosition,wxSize(200,20),wxTE_PROCESS_ENTER)
	{ this->s = s; this->i = i; }

	void OnValue(wxCommandEvent &event);

DECLARE_EVENT_TABLE()
};

class SurfaceTree : public wxSplitterWindow
{

public:
	Surface *selected;

	wxScrolledWindow *surflistpan;	// panel listing surfaces
	wxListCtrl *surflist;	// list of surfaces
	wxScrolledWindow *surfparpan;	// panel of parameters for current surface
	wxChoice *newsurf;	// which new surface to add


	SurfaceTree(wxWindow *parent, wxWindowID id = -1,
		wxPoint pos = wxDefaultPosition, wxSize size = wxDefaultSize);

	void DrawParamPanel(Surface *s);

	void Add(Surfaces *surfs);
	void Add(Surface *s);
	void OnRename(wxListEvent &event);
	void OnSelect(wxListEvent &event);
	void OnNew(wxCommandEvent &event);
	void OnKeyPress(wxListEvent &event);

	void OnUpDegree(wxCommandEvent &event) {
		wxButton *me = dynamic_cast<wxButton *>(event.GetEventObject());
		Algebraic *alg = (Algebraic *)me->GetClientData();
		alg->degree(alg->degree() + 1);
		surfparpan->DestroyChildren();
		DrawParamPanel(selected);
	}

	void OnDownDegree(wxCommandEvent &event) {
		wxButton *me = dynamic_cast<wxButton *>(event.GetEventObject());
		Algebraic *alg = (Algebraic *)me->GetClientData();
		alg->degree(alg->degree() - 1);
		surfparpan->DestroyChildren();
		DrawParamPanel(selected);
	}

DECLARE_EVENT_TABLE()
};

class SurfaceTreeParamButton : public wxButton
{
public:
	SurfaceTreeParamButton(wxWindow *parent, const wxString &text = "", const wxString &tooltip = "")
	: wxButton(parent,ID_SURFACETREEPARAMBUTTON,text) { SetToolTip(tooltip); cb = NULL; }
	
	SurfParamButton::Callback *cb;
	
	void OnSurfaceTreeParamButton(wxCommandEvent &event);
	
	DECLARE_EVENT_TABLE()
};

//Elmar:
//added combo boxes to surface params
class SurfaceTreeParamComboBox : public wxComboBox
{	
public:
	SurfaceTreeParamComboBox(	wxWindow * parent, const wxArrayString & choices, 
								const wxString & shortname="", 
								const wxString & desc="")
	: wxComboBox(	parent,ID_SURFACETREEPARAMCOMBOBOX, 
					choices[0], wxDefaultPosition, wxDefaultSize, 
					choices,wxCB_READONLY
				)
	  {SetToolTip(desc);cb=NULL;};
	SurfParamComboBox::Callback *cb;
	
	void OnSurfaceTreeComboBoxSelectingItem(wxCommandEvent &event);
	
	DECLARE_EVENT_TABLE()
};

#endif