#include "wxModeler.h"
#include "wx/colour.h"
#include "ParticleTwiddle.h"

BEGIN_EVENT_TABLE(ParticleTwiddle, wxScrolledWindow)
	EVT_BUTTON(ID_EXPAND,ParticleTwiddle::OnExpand)
	EVT_TEXT_ENTER(ID_RENAME,ParticleTwiddle::OnRename)
	EVT_SIZE(ParticleTwiddle::OnSize)
END_EVENT_TABLE()

void
ParticleTwiddle::OnExpand(wxCommandEvent &event)
{
	ParticleExpandButton *peb = (ParticleExpandButton *)(event.GetEventObject());
	peb->Toggle();
	topsizer->Layout();
}

void
ParticleTwiddle::OnRename(wxCommandEvent &event)
{
	ParticleName *pn = (ParticleName *)(event.GetClientObject());
	if (pn)
		pn->SetName(pn->GetValue());
}

void
ParticleTwiddle::OnSize(wxSizeEvent &event)
{
	topsizer->Layout();
}

ParticleTwiddle::ParticleTwiddle(wxWindow *parent)
: wxScrolledWindow(parent,-1)
{
	SetBackgroundColour(*wxWHITE);

	topsizer = new wxBoxSizer(wxVERTICAL);

	ParticleSystems *particlesystems = &(wxGetApp().particlesystems);
	for (int i = 0; i < particlesystems->size(); i++) {
		topsizer->Add(Add((*particlesystems)[i]));
	}

	SetSizerAndFit(topsizer);
	topsizer->Layout();
}

wxSizer *
ParticleTwiddle::Add(ParticleSystem *ps)
{
	wxBoxSizer *pssizer = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer *rowsizer = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *childsizer = new wxBoxSizer(wxVERTICAL);

	rowsizer->Add(new ParticleExpandButton(this,pssizer,childsizer),0,0,0);
	rowsizer->Add(new ParticleSystemName(this,ps),1,wxEXPAND | wxLEFT,4);

	pssizer->Add(rowsizer,0,0,0);
	pssizer->Add(childsizer,0,0,0);

	for (int i = 0; i < ps->particles.size(); i++) {
		childsizer->Add(Add(ps->particles[i]));
	}

	return pssizer;
}

wxSizer *
ParticleTwiddle::Add(Particles *p)
{
	wxBoxSizer *psizer = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer *rowsizer = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *childsizer = new wxBoxSizer(wxVERTICAL);

	rowsizer->Add(12,12);
	rowsizer->Add(new ParticleExpandButton(this,psizer,childsizer),0,0,0);
	rowsizer->Add(new ParticlesName(this,p),1,wxEXPAND | wxLEFT,4);

	psizer->Add(rowsizer,0,0,0);
	psizer->Add(childsizer,0,0,0);

	childsizer->Add(Add(p->attributes));
	childsizer->Add(Add(p->behaviors));
	childsizer->Add(Add(p->shaders));

	return psizer;
}

wxSizer *
ParticleTwiddle::Add(ParticleAttributes &attributes)
{
	wxBoxSizer *attrsizer = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer *rowsizer = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *childsizer = new wxBoxSizer(wxVERTICAL);

	rowsizer->Add(24,12);
	rowsizer->Add(new ParticleExpandButton(this,attrsizer,childsizer),0,0,0);
	rowsizer->Add(new wxTextCtrl(this,-1,_T("Attributes"),wxDefaultPosition,wxDefaultSize,wxTE_READONLY | wxNO_BORDER),1,wxEXPAND | wxLEFT,4);

	attrsizer->Add(rowsizer,0,0,0);
	attrsizer->Add(childsizer,0,0,0);

	for (ParticleAttributes::iterator attr = attributes.begin(); attr != attributes.end(); attr++) {
		childsizer->Add(Add((ParticleStuff *)(attr->second)));
	}

	return attrsizer;
}

wxSizer *
ParticleTwiddle::Add(ParticleBehaviors &behaviors)
{
	wxBoxSizer *behasizer = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer *rowsizer = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *childsizer = new wxBoxSizer(wxVERTICAL);

	rowsizer->Add(24,12);
	rowsizer->Add(new ParticleExpandButton(this,behasizer,childsizer),0,0,0);
	rowsizer->Add(new wxTextCtrl(this,-1,_T("Behaviors"),wxDefaultPosition,wxDefaultSize,wxTE_READONLY | wxNO_BORDER),1,wxEXPAND | wxLEFT,4);

	behasizer->Add(rowsizer,0,0,0);
	behasizer->Add(childsizer,0,0,0);

	for (ParticleBehaviors::iterator beha = behaviors.begin(); beha != behaviors.end(); beha++) {
		childsizer->Add(Add((ParticleStuff *)*beha));
	}

	return behasizer;
}

wxSizer *
ParticleTwiddle::Add(ParticleShaders &shaders)
{
	wxBoxSizer *shadsizer = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer *rowsizer = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *childsizer = new wxBoxSizer(wxVERTICAL);

	rowsizer->Add(24,12);
	rowsizer->Add(new ParticleExpandButton(this,shadsizer,childsizer),0,0,0);
	rowsizer->Add(new wxTextCtrl(this,-1,_T("Shaders"),wxDefaultPosition,wxDefaultSize,wxTE_READONLY | wxNO_BORDER),1,wxEXPAND | wxLEFT,4);

	shadsizer->Add(rowsizer,0,0,0);
	shadsizer->Add(childsizer,0,0,0);

	for (ParticleShaders::iterator shad = shaders.begin(); shad != shaders.end(); shad++) {
		childsizer->Add(Add((ParticleStuff *)*shad));
	}

	return shadsizer;
}

wxSizer *
ParticleTwiddle::Add(ParticleStuff *stuff)
{
	wxBoxSizer *stuffsizer = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer *rowsizer = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *childsizer = new wxBoxSizer(wxVERTICAL);

	rowsizer->Add(36,12);
	rowsizer->Add(new ParticleExpandButton(this,stuffsizer,childsizer),0,0,0);
	rowsizer->Add(new ParticleStuffName(this,stuff),1,wxEXPAND | wxLEFT,4);

	stuffsizer->Add(rowsizer,0,0,0);
	stuffsizer->Add(childsizer,0,0,0);

	return stuffsizer;
}