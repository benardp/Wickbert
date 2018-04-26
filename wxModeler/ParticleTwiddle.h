#ifndef PARTICLETWIDDLE_H
#define PARTICLETWIDDLE_H

#include "wx/wx.h"
#include "Particles/Particles.h"
#include "Particles/ParticleSystem.h"
#include "Particles/ParticleStuff.h"

enum {
	ID_EXPAND = 1,
	ID_RENAME,
};

class ParticleExpandButton : public wxButton
{
public:
	wxSizer *childsizer;
	wxSizer *parentsizer;
	bool childshown;

	ParticleExpandButton(wxWindow *parent, wxSizer *parentsizer, wxSizer *childsizer)
		:wxButton(parent,ID_EXPAND,"-",wxDefaultPosition,wxSize(12,12),wxNO_BORDER)
	{
		this->parentsizer = parentsizer;
		this->childsizer = childsizer;
		childshown = TRUE;
	}

	void Toggle() {
		childshown = !childshown;
		parentsizer->Show(childsizer,childshown);
		SetLabel(childshown ? _T("-") : _T("+"));
	}
};

class ParticleName : public wxTextCtrl
{
public:
	ParticleName(wxWindow *parent, wxString name)
		:wxTextCtrl(parent,ID_RENAME,name,wxDefaultPosition,wxDefaultSize,wxNO_BORDER) {}
	virtual void SetName(const wxString &name) {}
};

class ParticleSystemName : public ParticleName
{
public:
	ParticleSystem *ps;

	ParticleSystemName(wxWindow *parent, ParticleSystem *ps)
		:ParticleName(parent,ps->name.c_str())
	{
		this->ps = ps;
		SetClientData(ps);
	}
	void SetName(const wxString &name) { ps->name = name.c_str(); }
};

class ParticlesName : public ParticleName
{
public:
	Particles *p;

	ParticlesName(wxWindow *parent, Particles *p)
		:ParticleName(parent,p->name.c_str())
	{
		this->p = p;
	}
	void SetName(const wxString &name) { p->name = name.c_str(); }
};

class ParticleStuffName : public ParticleName
{
public:
	ParticleStuff *stuff;

	ParticleStuffName(wxWindow *parent, ParticleStuff *stuff)
		:ParticleName(parent,stuff->name.c_str())
	{
		this->stuff = stuff;
	}
	void SetName(const wxString &name) { stuff->name = name.c_str(); }
};

class ParticleTwiddle : public wxScrolledWindow
{
public:
	wxBoxSizer *topsizer;

	ParticleTwiddle(wxWindow *parent);
	wxSizer* Add(ParticleSystem *ps);
	wxSizer* Add(Particles *p);
	wxSizer *Add(ParticleAttributes &attr);
	wxSizer *Add(ParticleBehaviors &beha);
	wxSizer *Add(ParticleShaders &shad);
	wxSizer* Add(ParticleStuff *p);

	void OnRename(wxCommandEvent &event);
	void OnExpand(wxCommandEvent &event);
	void OnSize(wxSizeEvent &event);

DECLARE_EVENT_TABLE()
};

#endif