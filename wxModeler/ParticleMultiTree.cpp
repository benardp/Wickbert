#include "wxModeler.h"
#include "ParticleTree.h"

ParticleTree::ParticleTree(wxWindow *parent, wxWindowID id, wxPoint pos, wxSize size)
: wxTreeMultiCtrl(parent, id, pos, size)
{
	SetBackgroundColour(*wxWHITE);
	SetSpacingY(0);
	ParticleSystems *particlesystems = &(wxGetApp().particlesystems);
	Add(particlesystems);
}

void
ParticleTree::Add(ParticleSystems *ps)
{
	ParticleTreeItem root;

	for (int i = 0; i < ps->size(); i++) {
		root = AddRoot((*ps)[i]->name.c_str());
		Add(root,(*ps)[i]);
	}
}

void
ParticleTree::Add(ParticleTreeItem &parent, ParticleSystem *ps)
{
	ParticleTreeItem pnode;
	ParticleTreeItem titlenode;

	ParticleAttributes::iterator attr;
	ParticleBehaviors::iterator beha;
	ParticleShaders::iterator shad;

	Particles *p;

	for (int i = 0; i < ps->particles.size(); i++) {
		p = ps->particles[i];

		pnode = AppendNode(parent,p->name.c_str());

		titlenode = AppendNode(pnode,_T("Attributes"));
		for (attr = p->attributes.begin(); attr != p->attributes.end(); attr++) {
			Add(titlenode,(ParticleStuff *)(attr->second));
		}

		titlenode = AppendNode(pnode,_T("Behaviors"));
		for (beha = p->behaviors.begin(); beha != p->behaviors.end(); beha++) {
			Add(titlenode,(ParticleStuff *)*beha);
		}

		titlenode = AppendNode(pnode,_T("Shaders"));
		for (shad = p->shaders.begin(); shad != p->shaders.end(); shad++) {
			Add(titlenode,(ParticleStuff *)*shad);
		}
	}
}

void
ParticleTree::Add(ParticleTreeItem &parent, ParticleStuff *stuff)
{
	ParticleTreeItem stuffnode;

	stuffnode = AppendNode(parent,stuff->name.c_str());
	for (int i = 0; i < stuff->qlen(); i++) {
		AppendWindow(stuffnode,new ParticleTreeParam(this,stuff,i));
	}
}

BEGIN_EVENT_TABLE(ParticleTreeParam, wxPanel)
	EVT_TEXT_ENTER(ID_PARTICLEPARAM,ParticleTreeParam::OnValue)
END_EVENT_TABLE()

ParticleTreeParam::ParticleTreeParam(wxWindow *parent, ParticleStuff *stuff, int i)
: wxPanel(parent,-1)
{
	DoubleVector q;
	NameVector qnames;

	stuff->getq(q);
	stuff->qname(qnames);

	this->stuff = stuff;
	this->i = i;

	wxBoxSizer *rowsizer = new wxBoxSizer(wxHORIZONTAL);
	rowsizer->Add(new wxStaticText(this,-1,qnames[i].c_str(),wxDefaultPosition,wxSize(100,-1)));
	rowsizer->Add(
		new wxTextCtrl(this,ID_PARTICLEPARAM,wxString::Format("%g",q[i]),wxDefaultPosition,wxSize(100,-1),wxNO_BORDER),
		1,wxEXPAND
	);

	this->SetSizer(rowsizer);
}

void
ParticleTreeParam::OnValue(wxCommandEvent &event)
{
	DoubleVector q;

	wxTextCtrl *tc = (wxTextCtrl *)(event.GetEventObject());
	ParticleTreeParam *ptp = (ParticleTreeParam *)(tc->GetParent());

	ptp->stuff->getq(q);
	tc->GetValue().ToDouble(&(q[ptp->i]));
	ptp->stuff->setq(q);
}