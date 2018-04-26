#ifndef PARTICLETREE_H
#define PARTICLETREE_H

#include "wx/wx.h"
#include "wx/treectrl.h"
#include "wx/listctrl.h"
#include "Particles/Particles.h"
#include "Particles/ParticleSystem.h"
#include "Particles/ParticleStuff.h"
#include "Particles/ParticleStuffParameters.h"

typedef wxTreeItemId ParticleTreeItem;

class ParticleTreeItemData : public wxTreeItemData
{
public:
	bool editable;
	bool draggable;

	virtual wxString GetName() = 0;
	virtual wxString SetName(const wxString &name) = 0;
	virtual bool Remove() { return false; }
	virtual void MakePanel(wxScrolledWindow *panel) {}
	ParticleTreeItemData() { editable = false; draggable = false; }
	virtual bool dropped(ParticleTreeItemData *dragged) { return false; }
};

class ParticleTreeHeader : public ParticleTreeItemData
{
public:
	virtual wxString GetName() { return _T(""); }
	virtual wxString SetName(const wxString &name) { return name; }
};

class ParticleTreeStuffHeader : public ParticleTreeItemData
{
public:
	wxString prefix;
	Particles *p;

	ParticleTreeStuffHeader(const wxString &prefix, Particles *p) {
		this->prefix = prefix; this->p = p;
	}
	virtual wxString GetName() { return _T(""); }
	virtual wxString SetName(const wxString &name) { return name; }
	virtual void MakePanel(wxScrolledWindow *panel);
};

class ParticleSystemItemData : public ParticleTreeItemData
{
public:
	ParticleSystem *ps;

	ParticleSystemItemData(ParticleSystem *ps) { editable = true; this->ps = ps; }
	wxString GetName() { return ps->name.c_str(); }
	wxString SetName(const wxString &name) { 
		ps->name = name.c_str();
		ps->attachAttributes();	// in case something references the name we just changed
		return name;
	}
	virtual void MakePanel(wxScrolledWindow *panel);
	bool Remove();
};

class ParticleSystemsItemData : public ParticleTreeItemData
{
public:
	ParticleSystemsItemData() { }
	wxString GetName() { return ""; }
	wxString SetName(const wxString &name) { return ""; }
	virtual void MakePanel(wxScrolledWindow *panel);
};

class ParticlesReport : public wxStaticText
{
public:
	Particles *p;

	ParticlesReport(wxWindow *parent, Particles *p)
	: wxStaticText(parent, ID_PARTICLESREPORT, "???")
	{ this->p = p; }

	void OnUpdateUI(wxUpdateUIEvent &event)
	{ SetLabel(wxString::Format("%d",p->size())); }

	DECLARE_EVENT_TABLE()
};

class ParticlesItemData : public ParticleTreeItemData
{
public:
	Particles *p;

	ParticlesItemData(Particles *p) { editable = true; this->p = p; }
	wxString GetName() { return p->name.c_str(); }
	wxString SetName(const wxString &name) {
		p->name = name.c_str();
		p->particleSystem->attachAttributes(); // in case something references the name we just changed
		return name;
	}
	bool Remove() {	return p->particleSystem->findanderase(p) == -1; }
	virtual void MakePanel(wxScrolledWindow *panel);
};

class ParticleStuffItemData : public ParticleTreeItemData
{
public:
	ParticleStuff *stuff;

	ParticleStuffItemData(ParticleStuff *stuff) {
		this->stuff = stuff; editable = true; draggable = true;
	}
	wxString GetName() { return stuff->name.c_str(); }
	wxString SetName(const wxString &name);
	bool Remove() {stuff->ps->remove(stuff); return true;}
	virtual void MakePanel(wxScrolledWindow *panel);
	virtual bool dropped(ParticleTreeItemData *dragged);
};

class ImplicitInterrogatorItemData : public ParticleStuffItemData
{
public:
	ImplicitInterrogatorItemData(ParticleStuff *stuff) : ParticleStuffItemData(stuff) { }
	virtual void MakePanel(wxScrolledWindow *panel);
};

class ParticleTreeParam : public wxTextCtrl
{
public:
	ParticleTreeParam(wxWindow *parent, const wxString &text = "")
	: wxTextCtrl(parent,ID_PARTICLETREEPARAM,text,wxDefaultPosition,wxSize(200,20),wxTE_PROCESS_ENTER) { }

	virtual void SetParam(const wxString &text) { }
	void OnParamChange(wxCommandEvent &event);
	void OnParamFocus(wxFocusEvent &event);

DECLARE_EVENT_TABLE()
};

class ParticleTreeParamButton : public wxButton
{
public:
	ParticleTreeParamButton(wxWindow *parent, const wxString &text = "", const wxString &tooltip = "")
	: wxButton(parent,ID_PARTICLETREEPARAMBUTTON,text) { SetToolTip(tooltip); cb = NULL; }

	PSParamButton::Callback *cb;

	void OnParticleTreeParamButton(wxCommandEvent &event);

DECLARE_EVENT_TABLE()
};

class ParticleTreeParamCheckBox : public wxCheckBox
{
public:
	PSParam *param;

	ParticleTreeParamCheckBox(wxWindow *parent, PSParam *param)
	: wxCheckBox(parent,ID_PARTICLETREEPARAMCHECKBOX,"(" + param->shortname + ")") {
		this->param = param;
		SetValue(*((bool *)param->ref()));
	}

	void OnUpdateUI(wxUpdateUIEvent &event) {
		SetValue(*((bool *)param->ref()));
	}

	void OnParticleTreeParamCheckBox(wxCommandEvent &event) {
		param->set(GetValue() ? "true" : "false");
	}

	DECLARE_EVENT_TABLE()
};

class ParticleTreeStuffParam : public ParticleTreeParam
{
public:
	ParticleStuff *stuff;
	int i;
	PSParam *param;

	ParticleTreeStuffParam(wxWindow *parent, const wxString &text, ParticleStuff *stuff, int i)
		: ParticleTreeParam(parent, text)
	{ this->stuff = stuff; this->i = i; }

	ParticleTreeStuffParam(wxWindow *parent, PSParam *param)
		: ParticleTreeParam(parent)
	{
		this->param = param;
		SetValue(param->get());
	}

	void SetParam(const wxString &text) {
		param->set(text.ToStdString());
	}

	void OnUpdateUI(wxUpdateUIEvent &event) {
		if (FindFocus() != this) SetValue(param->get());
	}

	DECLARE_EVENT_TABLE()
};

class ParticleTreeStuffParamPerParticle : public ParticleTreeParam
{
public:
	ParticleStuff *stuff;
	int i;
	PSParamPerParticle *param;
	int selected;

	ParticleTreeStuffParamPerParticle(wxWindow *parent, const wxString &text, ParticleStuff *stuff, int i, int selected)
		: ParticleTreeParam(parent, text)
	{ this->stuff = stuff; this->i = i; this->selected = selected; }

	ParticleTreeStuffParamPerParticle(wxWindow *parent, PSParamPerParticle *param, int selected)
		: ParticleTreeParam(parent)
	{
		this->param = param;
		this->selected = selected;
		SetValue(selected == -1 ? "<unselected>" : param->get(selected).c_str());
	}

	void SetParam(const wxString &text) {
		if (selected >= 0 && selected < (int) param->stuff->ps->size())
			param->set(selected,text.ToStdString());
	}

	void OnUpdateUI(wxUpdateUIEvent &event) {
		// Don't update if currently entering data
		if (FindFocus() != this) {
			selected = param->stuff->ps->selectedParticle;
			SetValue(selected == -1 ? "<unselected>" : param->get(selected).c_str());
		}
	}

DECLARE_EVENT_TABLE()
};

class ParticleTreeShowSelected : public ParticleTreeParam
{
public:
	Particles *p;

	ParticleTreeShowSelected(wxWindow *parent, Particles *p)
		: ParticleTreeParam(parent)
	{ this->p = p; UpdateNumber(); }

	void OnUpdateUI(wxUpdateUIEvent &event) { UpdateNumber(); }

	void SetParam(const wxString &text) {
		p->selectedParticle = atoi(text.c_str());
		if (p->selectedParticle < 0 || p->selectedParticle >= (int) p->size()) {
			p->selectedParticle = -1;
			SetValue("<unselected>");
		}
	}
		
	void UpdateNumber() {
		if (FindFocus() != this)
			SetValue(p->selectedParticle == -1 ? "<unselected>" : wxString::Format("%d",p->selectedParticle));
	}

DECLARE_EVENT_TABLE()
};

class ParticleTreeStuffAA : public ParticleTreeParam
{
public:
	ParticleStuff *stuff;
	AttachedAttribute *aa;

	ParticleTreeStuffAA(wxWindow *parent, AttachedAttribute *aa)
		: ParticleTreeParam(parent)
	{
		this->aa = aa;
		SetValue(aa->attr_name.c_str());
	}

	void SetParam(const wxString &text) {
		aa->attr_name = text.c_str();
		aa->attach();
	}
};

class ParticleTree : public wxTreeCtrl
{
public:
	wxScrolledWindow *parampanel;

	std::map<ParticleStuff *, int> iconmap;
	wxImageList *imagelist;

	ParticleTree(wxWindow *parent, wxWindowID id = -1,
		wxPoint pos = wxDefaultPosition, wxSize size = wxDefaultSize);

	void Reload();

	void Add(ParticleSystems *ps);
	void Add(ParticleTreeItem &parent, ParticleSystem *ps);
	void Add(ParticleTreeItem &parent, Particles *p);
	void Add(ParticleTreeItem &parent, const ParticleAttributes &attrs);
	void Add(ParticleTreeItem &parent, const ParticleBehaviors &behas);
	void Add(ParticleTreeItem &parent, const ParticleShaders &shads);
	void Add(ParticleTreeItem &parent, ParticleStuff *stuff);

	void OnEdit(wxTreeEvent &event);
	void OnRename(wxTreeEvent &event);
	void OnTreeSelect(wxTreeEvent &event);
	void OnDelete(wxTreeEvent &event);
	void OnKeyPress(wxTreeEvent &event);
	// void OnParamChange(wxCommandEvent &event);

	ParticleTreeItemData *dragged;

	void OnBeginDrag(wxTreeEvent &event);
	void OnEndDrag(wxTreeEvent &event);

DECLARE_EVENT_TABLE()
};

class ParticleWidget : public wxSplitterWindow
{
public:
	ParticleTree *partree;
	wxScrolledWindow *parpanel;

	ParticleWidget(wxWindow *parent, wxWindowID id = -1,
		wxPoint pos = wxDefaultPosition, wxSize size = wxDefaultSize)
		: wxSplitterWindow(parent,-1,wxDefaultPosition,wxDefaultSize,wxSP_3D)
	{
		Reload();
	}

	void Reload() {
		DestroyChildren();
		partree = new ParticleTree(this,ID_PARTICLETREE);
		parpanel = new wxScrolledWindow(this,-1,wxDefaultPosition,wxDefaultSize,wxSIMPLE_BORDER);
		wxSizer *sizer = new wxBoxSizer(wxVERTICAL);
		parpanel->SetSizer(sizer);
		partree->parampanel = parpanel;
		//split this after the vertical split:
		//SplitHorizontally(partree,parpanel); 
	}

	void OnImpIntChanged(wxCommandEvent &event);
	void OnNewParticleSystem(wxCommandEvent &event);
	void OnNewParticles(wxCommandEvent &event);
	void OnNewParticleStuff(wxCommandEvent &event);
	void OnCopyParticles(wxCommandEvent &event);
	void OnAddParticle(wxCommandEvent &event);
	void OnKillParticles(wxCommandEvent &event);
	void OnPlayParticles(wxCommandEvent &event);
	void OnShowParticles(wxCommandEvent &event);

DECLARE_EVENT_TABLE()
};

class ParticleTreeParamComboBox : public wxComboBox
{	
public:
	ParticleTreeParamComboBox(	wxWindow * parent, const wxArrayString & choices, 
								const wxString & shortname="", 
								const wxString & desc="")
	: wxComboBox(	parent,ID_PARTICLEPARAMCOMBOBOX, 
					choices[0], wxDefaultPosition, wxDefaultSize, 
					choices,wxCB_READONLY
				)
	  {SetToolTip(desc);cb=NULL;};

	PSParamComboBox::Callback *cb;
	
	void OnParticleTreeComboBoxSelectingItem(wxCommandEvent &event);
	
	DECLARE_EVENT_TABLE()
};

#endif

