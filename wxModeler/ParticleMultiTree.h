#ifndef PARTICLETREE_H
#define PARTICLETREE_H

#include "wx/wx.h"
#include "wx/wxTreeMultiCtrl.h"
#include "Particles/Particles.h"
#include "Particles/ParticleSystem.h"
#include "Particles/ParticleStuff.h"

#if 0
class ParticleTreeItem : public wxTreeMultiItem
{
};
#else
typedef wxTreeMultiItem ParticleTreeItem;
#endif

class ParticleTreeParam : public wxPanel
{
public:
	ParticleStuff *stuff;
	int i;

	ParticleTreeParam(wxWindow *parent, ParticleStuff *stuff, int i);
	void OnValue(wxCommandEvent &event);

DECLARE_EVENT_TABLE()
};

class ParticleTree : public wxTreeMultiCtrl
{
public:
	ParticleTree(wxWindow *parent, wxWindowID id = -1,
		wxPoint pos = wxDefaultPosition, wxSize size = wxDefaultSize);

	void Add(ParticleSystems *systems);
	void Add(ParticleTreeItem &parent, ParticleSystem *ps);
	void Add(ParticleTreeItem &parent, ParticleStuff *stuff);
};

#endif