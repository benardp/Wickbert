#include "wxModeler.h"
#include "ParticleTree.h"
#include "wx/listctrl.h"
#include "wx/imaglist.h"
#include "Particles/Attributes/ImplicitInterrogator.h"
#include "Particles/ParticleShader.h"
#include "Particles/ParticleBehavior.h"

#include "DisplayWindow.h" // for the a-z zoom hack

BEGIN_EVENT_TABLE(ParticleTree, wxTreeCtrl)
	EVT_TREE_BEGIN_LABEL_EDIT(ID_PARTICLETREE,ParticleTree::OnEdit)
	EVT_TREE_END_LABEL_EDIT(ID_PARTICLETREE,ParticleTree::OnRename)
	EVT_TREE_SEL_CHANGED(ID_PARTICLETREE,ParticleTree::OnTreeSelect)
	EVT_TREE_KEY_DOWN(ID_PARTICLETREE,ParticleTree::OnKeyPress)
	EVT_TREE_BEGIN_DRAG(ID_PARTICLETREE, ParticleTree::OnBeginDrag)
	EVT_TREE_END_DRAG(ID_PARTICLETREE, ParticleTree::OnEndDrag)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(ParticleWidget, wxSplitterWindow)
	EVT_CHOICE(ID_IMPINTCHANGED,ParticleWidget::OnImpIntChanged)
	EVT_CHOICE(ID_NEWPARTICLESTUFF,ParticleWidget::OnNewParticleStuff)
	EVT_CHOICE(ID_COPYPARTICLES,ParticleWidget::OnCopyParticles)
	EVT_TEXT_ENTER(ID_NEWPARTICLES,ParticleWidget::OnNewParticles)
	EVT_TEXT_ENTER(ID_NEWPARTICLESYSTEM,ParticleWidget::OnNewParticleSystem)
	EVT_BUTTON(ID_ADDPARTICLE,ParticleWidget::OnAddParticle)
	EVT_BUTTON(ID_KILLPARTICLES,ParticleWidget::OnKillParticles)
	EVT_BUTTON(ID_PLAYPARTICLES,ParticleWidget::OnPlayParticles)
	EVT_BUTTON(ID_SHOWPARTICLES,ParticleWidget::OnShowParticles)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(ParticleTreeParam, wxTextCtrl)
	EVT_TEXT_ENTER(ID_PARTICLETREEPARAM,ParticleTreeParam::OnParamChange)
	EVT_KILL_FOCUS(ParticleTreeParam::OnParamFocus)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(ParticleTreeParamButton, wxButton)
	EVT_BUTTON(ID_PARTICLETREEPARAMBUTTON,ParticleTreeParamButton::OnParticleTreeParamButton)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(ParticleTreeParamComboBox, wxComboBox)
EVT_COMBOBOX(ID_PARTICLEPARAMCOMBOBOX, ParticleTreeParamComboBox::OnParticleTreeComboBoxSelectingItem)
END_EVENT_TABLE()


BEGIN_EVENT_TABLE(ParticleTreeParamCheckBox, wxCheckBox)
	EVT_CHECKBOX(ID_PARTICLETREEPARAMCHECKBOX,ParticleTreeParamCheckBox::OnParticleTreeParamCheckBox)
	EVT_UPDATE_UI(wxID_ANY,ParticleTreeParamCheckBox::OnUpdateUI)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(ParticleTreeStuffParam, ParticleTreeParam)
	EVT_UPDATE_UI(wxID_ANY,ParticleTreeStuffParam::OnUpdateUI)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(ParticleTreeStuffParamPerParticle, ParticleTreeParam)
	EVT_UPDATE_UI(wxID_ANY,ParticleTreeStuffParamPerParticle::OnUpdateUI)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(ParticleTreeShowSelected, ParticleTreeParam)
	EVT_UPDATE_UI(wxID_ANY,ParticleTreeShowSelected::OnUpdateUI)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(ParticlesReport, wxStaticText)
	EVT_UPDATE_UI(ID_PARTICLESREPORT,ParticlesReport::OnUpdateUI)
END_EVENT_TABLE()


void
ParticleTree::OnBeginDrag(wxTreeEvent &event)
{
	wxTreeItemId id = event.GetItem();
	dragged = (ParticleTreeItemData *)(GetItemData(id));
	if (dragged && dragged->draggable)
		event.Allow();
}

void
ParticleTree::OnEndDrag(wxTreeEvent &event)
{
// I think this has changed in 2.6.0 -jch

#if 1
	wxTreeItemId dropped_id = event.GetItem();
	ParticleTreeItemData *dropped = (ParticleTreeItemData *)(GetItemData(dropped_id));
	if (dropped->dropped(dragged)) {
		wxTreeItemId dragged_id = dragged->GetId();
		ParticleStuffItemData *psid = (ParticleStuffItemData *)(dragged);
		InsertItem(GetItemParent(dropped_id),dropped_id,GetItemText(dragged_id),GetItemImage(dragged_id),-1,new ParticleStuffItemData(psid->stuff));
		Delete(dragged_id);
	}
#endif
}

/** Assumed dragged is a ParticleStuffItemData and has a stuff pointer.
 * 
 */
bool ParticleStuffItemData::dropped(ParticleTreeItemData *dragged)
{
#if 0
	ParticleStuffItemData *psid = dynamic_cast<ParticleStuffItemData *>(dragged);
	if (!psid) return false;
#else
	ParticleStuffItemData *psid = (ParticleStuffItemData *)(dragged);
#endif

	return psid->stuff->moveAfter(stuff);
}

void
ParticleTree::OnEdit(wxTreeEvent &event)
{
	wxTreeItemId id = event.GetItem();
	ParticleTreeItemData *ptid = (ParticleTreeItemData *)(GetItemData(id));
	wxString newname = event.GetLabel();
	if (!ptid->editable) {
		event.Veto();
		return;
	}

	/* Remove the [classname] added by SetName when the name != classname */
	if (newname.Last() == ']')
		SetItemText(id,newname.BeforeLast(' '));
}

/* I really wanted this to change the name to "name [classname]" after editing
 * but no matter what I tried, I couldn't change the name assigned to the item
 * in this callback.
 */
void
ParticleTree::OnRename(wxTreeEvent &event)
{
	wxString newname;

	wxTreeItemId id = event.GetItem();
	ParticleTreeItemData *ptid = (ParticleTreeItemData *)(GetItemData(id));

	/* if item is cancelled, still need to add back the [classname]
	 * that we removed in the OnEdit routine */
	if (event.IsEditCancelled())
		newname = GetItemText(id);
	else
		newname = event.GetLabel();
	SetItemText(id,ptid->SetName(newname));

	/* If we don't veto, then the user-changed label overwrites the label
	   we just wrote to the ParticleTree */
	event.Veto();
}

void
ParticleTree::OnKeyPress(wxTreeEvent &event)
{
	switch (event.GetKeyCode()) {
		case WXK_DELETE:
		//added for wxWidgets 2.6.3 (I don't know why 
		//Delete key generates event with key code
		// 389 instead of WXK_DELETE) ---ms
		case 389:
			{
			wxTreeItemId id = GetSelection(); // event.GetItem() doesn't seem to work here
			ParticleTreeItemData *ptid = (ParticleTreeItemData *)(GetItemData(id));

			ParticleSystems *particlesystems = &(wxGetApp().particlesystems);
			ParticleSystems::iterator iter;
			std::vector<Particles*>::iterator piter;
			std::map<std::string, ParticleAttribute*>::iterator iterAttr;
			std::vector<ParticleBehavior*>::iterator iterBehav;
			std::vector<ParticleShader*>::iterator iterShad;
			std::vector<PSParam *>::iterator iterParam;
			AttachedAttributes::iterator iterAtAt;
			//iterate over all ParticleStuff objects of all Particles objects to set to null the deleted ParticleStuff
			if (ptid->Remove()) {
				//iterate over particle systems
				for(iter = particlesystems->begin(); iter != particlesystems->end(); iter++){
					//iterate over particle collections
					for( piter = ((*iter)->particles).begin();  piter != ((*iter)->particles).end(); piter++){
						//iter over stuff
						for( iterAttr = ((*piter)->attributes).begin();  iterAttr != ((*piter)->attributes).end(); iterAttr++){
							//iter over params
							for( iterParam = (iterAttr->second->params).begin(); iterParam != (iterAttr->second->params).end(); iterParam++){
								(*iterParam)->set((*iterParam)->get());
							}
							for(iterAtAt = iterAttr->second->attachedattributes.begin(); iterAtAt != iterAttr->second->attachedattributes.end(); iterAtAt++){
								(*iterAtAt)->attach();
							}
						}
						for( iterBehav = ((*piter)->behaviors).begin();  iterBehav != ((*piter)->behaviors).end(); iterBehav++){
							//iter over params
							for( iterParam = ((*iterBehav)->params).begin(); iterParam != ((*iterBehav)->params).end(); iterParam++){
								(*iterParam)->set((*iterParam)->get());
							}
							for(iterAtAt = (*iterBehav)->attachedattributes.begin(); iterAtAt != (*iterBehav)->attachedattributes.end(); iterAtAt++){
								(*iterAtAt)->attach();
							}
						}
						for( iterShad = ((*piter)->shaders).begin();  iterShad != ((*piter)->shaders).end(); iterShad++){
							//iter over params
							for( iterParam = ((*iterShad)->params).begin(); iterParam != ((*iterShad)->params).end(); iterParam++){
								(*iterParam)->set((*iterParam)->get());
							}
							for(iterAtAt = (*iterShad)->attachedattributes.begin(); iterAtAt != (*iterShad)->attachedattributes.end(); iterAtAt++){
								(*iterAtAt)->attach();
							}
						}
					}
				}
			}	
			Delete(id); /* Delete() is a TreeCtrl call */
			particlesystems->attachAttributes();
			particlesystems->surfacesAttach();
			}
			break;
		case WXK_F5:
			Reload();
			break;
	}
}

bool
ParticleSystemItemData::Remove()
{
	ParticleSystems *particlesystems = &(wxGetApp().particlesystems);
	ParticleSystems::iterator psit = find(particlesystems->begin(), particlesystems->end(), ps);
	if (psit == particlesystems->end()) return false;
	particlesystems->erase(psit);
	return true;
}

void
ParticleTreeParam::OnParamChange(wxCommandEvent &event)
{
	ParticleTreeParam *ptp = (ParticleTreeParam *)(event.GetEventObject());
	wxString val = event.GetString();
	ptp->SetParam(val);
	//[tk 7.12.05] Particles reference parameters attach() should go here.  Where one attaches, all should re-attach
}

void
ParticleTreeParam::OnParamFocus(wxFocusEvent &event)
{
	ParticleTreeParam *ptp = (ParticleTreeParam *)(event.GetEventObject());
	wxString val = ptp->GetValue();
	ptp->SetParam(val);
}

#include "bitmaps/ParticleSystems.xpm"
#include "bitmaps/ParticleSystem.xpm"
#include "bitmaps/Particles.xpm"

static char *stuffheader_xpm[] = {
	"16 16 3 1", /* columns, rows, colors, chars-per-pixel */
	"  c None",
	"r c #FF0000",
	"k c #000000",
	/* pixels */
	"                ",
	"                ",
	"        kkkk    ",
	" kkkkkkk    kk  ",
	"k       kkkk  k ",
	"k             k ",
	"k             k ",
	"k             k ",
	"k             k ",
	"k             k ",
	"k             k ",
	"k             k ",
	"k             k ",
	" kkkkkkkkkkkkk  ",
	"                ",
	"                "
};

static char *stuff_xpm[] = {
	"16 16 3 1", /* columns, rows, colors, chars-per-pixel */
	"  c None",
	"r c #FF0000",
	"k c #000000",
	/* pixels */
	"                ",
	"                ",
	" kkkkkkkkkkkkk  ",
	" k           k  ",
	" k  kkkkkkk  k  ",
	" k           k  ",
	" k  kkkkkkk  k  ",
	" k           k  ",
	" k  kkkkkkk  k  ",
	" k           k  ",
	" k  kkkkkkk  k  ",
	" k           k  ",
	" k           k  ",
	" kkkkkkkkkkkkk  ",
	"                ",
	"                "
};

ParticleTree::ParticleTree(wxWindow *parent, wxWindowID id, wxPoint pos, wxSize size)
: wxTreeCtrl(parent, id, pos, size, wxTR_HAS_BUTTONS | wxTR_EDIT_LABELS)
{
#if 0
	imagelist = new wxImageList(16,15);

	wxBitmap *psystems_bitmap = new wxBitmap("bitmaps/ParticleSystems.bmp",wxBITMAP_TYPE_BMP);
	wxBitmap *ps_bitmap = new wxBitmap("bitmaps/ParticleSystem.bmp",wxBITMAP_TYPE_BMP);
	wxBitmap *p_bitmap = new wxBitmap("bitmaps/Particles.bmp",wxBITMAP_TYPE_BMP);
	wxBitmap *stuffheader_bitmap = new wxBitmap("bitmaps/ParticleStuffHeader.bmp",wxBITMAP_TYPE_BMP);
	wxBitmap *stuff_bitmap = new wxBitmap("bitmaps/ParticleStuff.bmp",wxBITMAP_TYPE_BMP);

	imagelist->Add(*psystems_bitmap);
	imagelist->Add(*ps_bitmap);
	imagelist->Add(*p_bitmap);
	imagelist->Add(*stuffheader_bitmap);
	imagelist->Add(*stuff_bitmap);

	AssignImageList(imagelist);
#else
	wxIcon icons[5];
    icons[0] = wxIcon(ParticleSystems_xpm);
    icons[1] = wxIcon(ParticleSystem_xpm);
    icons[2] = wxIcon(Particles_xpm);
    icons[3] = wxIcon(stuffheader_xpm);
    icons[4] = wxIcon(stuff_xpm);

	imagelist = new wxImageList(icons[0].GetWidth(),icons[0].GetHeight());

	imagelist->Add(icons[0]);
	imagelist->Add(icons[1]);
	imagelist->Add(icons[2]);
	imagelist->Add(icons[3]);
	imagelist->Add(icons[4]);

	AssignImageList(imagelist);
#endif
	SetBackgroundColour(*wxWHITE);
	Reload();
}

void
ParticleTree::Reload()
{
	DeleteAllItems();
	ParticleSystems *particlesystems = &(wxGetApp().particlesystems);
	ParticleTreeItem root = AddRoot(_T("Particle Systems"),0,-1,new ParticleSystemsItemData);
	Add(particlesystems);
}

void
ParticleTree::Add(ParticleSystems *ps)
{
	ParticleTreeItem root = GetRootItem();

	for (unsigned int i = 0; i < ps->size(); i++) {
		Add(root,(*ps)[i]);
	}

	Expand(root);
}

void
ParticleTree::Add(ParticleTreeItem &parent, ParticleSystem *ps)
{
	ParticleTreeItem psnode = AppendItem(parent,ps->name.c_str(), 1, -1, new ParticleSystemItemData(ps));

	for (unsigned int i = 0; i < ps->particles.size(); i++)
		Add(psnode,ps->particles[i]);

	Expand(psnode);
}

void
ParticleTree::Add(ParticleTreeItem &parent, Particles *p)
{
	ParticleTreeItem pnode;
	ParticleTreeItem titlenode;

	pnode = AppendItem(parent,p->name.c_str(), 2, -1, new ParticlesItemData(p));

	titlenode = AppendItem(pnode,_T("Attributes"),3,-1,new ParticleTreeStuffHeader(_T("Attribute:"),p));
	Add(titlenode,p->attributes);

	titlenode = AppendItem(pnode,_T("Behaviors"),3,-1,new ParticleTreeStuffHeader(_T("Behavior:"),p));
	Add(titlenode,p->behaviors);

	titlenode = AppendItem(pnode,_T("Shaders"),3,-1,new ParticleTreeStuffHeader(_T("Shader:"),p));
	Add(titlenode,p->shaders);

	Expand(pnode);
}

void
ParticleTree::Add(ParticleTreeItem &parent, const ParticleAttributes &attrs)
{
	ParticleAttributes::const_iterator attr;
	for (attr = attrs.begin(); attr != attrs.end(); attr++) {
		Add(parent,(ParticleStuff *)(attr->second));
	}
	Expand(parent);
}

void
ParticleTree::Add(ParticleTreeItem &parent, const ParticleBehaviors &behas)
{
	ParticleBehaviors::const_iterator beha;
	for (beha = behas.begin(); beha != behas.end(); beha++) {
		Add(parent,(ParticleStuff *)*beha);
	}
	Expand(parent);
}

void
ParticleTree::Add(ParticleTreeItem &parent, const ParticleShaders &shads)
{
	ParticleShaders::const_iterator shad;
	for (shad = shads.begin(); shad != shads.end(); shad++) {
		Add(parent,(ParticleStuff *)*shad);
	}
	Expand(parent);
}

void
ParticleTree::Add(ParticleTreeItem &parent, ParticleStuff *stuff)
{
	ParticleTreeItem stuffnode;
	wxString stuffname;

	if (wxString(stuff->getName().c_str()) == wxString(stuff->getClass().c_str()).AfterLast(':')) {
		stuffname = stuff->getName().c_str();
	} else {
		stuffname.Printf("%s [%s]",stuff->getName().c_str(),stuff->getClass().c_str());
	}

	ParticleStuffItemData *psid;
	if (false && dynamic_cast<ImplicitInterrogator *>(stuff)) {
		psid = new ImplicitInterrogatorItemData(stuff);
	} else {
		psid = new ParticleStuffItemData(stuff);
	}

#if 0
	wxBitmap *bitmap = new wxBitmap(stuff->xpm);
	iconmap[stuff] = imagelist->Add(*bitmap);
#endif

	stuffnode = AppendItem(parent, stuffname, 4, -1, psid);
}

wxString
ParticleStuffItemData::SetName(const wxString &name)
{
	wxString stuffname = wxString(stuff->getClass().c_str()).AfterLast(':');

	stuff->setName(name.ToStdString());

	if (stuff->getName().c_str() != stuffname) {
//		stuffname.Printf("%s [%s]",stuff->getName().c_str(),stuffname);
		stuffname = wxString(stuff->getName().c_str()) + _T(" [") + stuffname + _T("]");
	}

	stuff->ps->attachAttributes();	// in case something referenced this name

	return stuffname;
}

void
ParticleTree::OnTreeSelect(wxTreeEvent &event)
{
	parampanel->DestroyChildren();

//	wxTreeItemId id = event.GetItem();
	wxTreeItemId id = GetSelection();
	ParticleTreeItemData *ptid = (ParticleTreeItemData *)(GetItemData(id));
	if (ptid != NULL)
		ptid->MakePanel(parampanel);

	Refresh();
}

void
ParticleSystemsItemData::MakePanel(wxScrolledWindow *panel)
{
	// make a textbox to create a new particle system

	wxSizer *rowsizer = new wxBoxSizer(wxHORIZONTAL);
	rowsizer->Add(new wxStaticText(panel,-1,_T("New ")),0,0,0);
	rowsizer->Add(new wxTextCtrl(panel,ID_NEWPARTICLESYSTEM, "", wxDefaultPosition, wxDefaultSize, wxTE_PROCESS_ENTER),1,0,0);

	wxSizer *sizer = new wxBoxSizer(wxVERTICAL);
	sizer->Add(rowsizer,0,wxEXPAND,0);

	panel->SetSizer(sizer);
	panel->Layout();
}

void
ParticleSystemItemData::MakePanel(wxScrolledWindow *panel)
{
	// make a textbox to create a new particles object

	wxSizer *rowsizer = new wxBoxSizer(wxHORIZONTAL);
	rowsizer->Add(new wxStaticText(panel,-1,_T("New ")),0,0,0);
	rowsizer->Add(new wxTextCtrl(panel,ID_NEWPARTICLES, "", wxDefaultPosition, wxDefaultSize, wxTE_PROCESS_ENTER),1,0,0);

	wxSizer *sizer = new wxBoxSizer(wxVERTICAL);
	sizer->Add(rowsizer,0,wxEXPAND,0);

	panel->SetSizer(sizer);
	panel->Layout();
}

void
ParticleWidget::OnAddParticle(wxCommandEvent &event)
{
	wxButton *addbutton = dynamic_cast<wxButton *>(event.GetEventObject());
	if (!addbutton) return;
	Particles *p = (Particles *)addbutton->GetClientData();
	p->addParticle();
}

void
ParticleWidget::OnKillParticles(wxCommandEvent &event)
{
	wxButton *addbutton = dynamic_cast<wxButton *>(event.GetEventObject());
	if (!addbutton) return;
	Particles *p = (Particles *)addbutton->GetClientData();
	p->removeAll();
}

void
ParticlesItemData::MakePanel(wxScrolledWindow *panel)
{
	wxFlexGridSizer *flexsizer = new wxFlexGridSizer(0,2,1,10);
	flexsizer->AddGrowableCol(0);
	flexsizer->AddGrowableCol(1);

	// Report number of particles
	flexsizer->Add(new wxStaticText(panel,-1,"# Particles:",wxDefaultPosition,wxDefaultSize,wxALIGN_RIGHT),1,wxEXPAND);
	flexsizer->Add(new ParticlesReport(panel,p),1,wxEXPAND);

	// Add particle button
	flexsizer->Add(new wxStaticText(panel,-1,"Add Particle:",wxDefaultPosition,wxDefaultSize,wxALIGN_RIGHT),1,wxEXPAND);
	wxButton *addbutton = new wxButton(panel,ID_ADDPARTICLE,"Add");
	addbutton->SetClientData(p);
	flexsizer->Add(addbutton,1,wxEXPAND);

	// Kill Particles button
	flexsizer->Add(new wxStaticText(panel,-1,"Kill Particles:",wxDefaultPosition,wxDefaultSize,wxALIGN_RIGHT),1,wxEXPAND);
	wxButton *killbutton = new wxButton(panel,ID_KILLPARTICLES,"Kill");
	killbutton->SetClientData(p);
	flexsizer->Add(killbutton,1,wxEXPAND);

	// make a textbox to copy attributes, behaviors and shaders from another particles
	flexsizer->Add(new wxStaticText(panel,-1,_T("Copy From:"),wxDefaultPosition,wxDefaultSize,wxALIGN_RIGHT),1,wxEXPAND);
	wxChoice *copyfrom = new wxChoice(panel, ID_COPYPARTICLES);
	for (unsigned int i = 0; i < p->particleSystem->particles.size(); i++) {
		Particles *pi = p->particleSystem->particles[i];
		if (pi != p)
			copyfrom->Append(pi->name.c_str(),pi);
	}
	flexsizer->Add(copyfrom,1,wxEXPAND);

	panel->SetSizer(flexsizer);
	panel->Layout();
}

void
ParticleTreeStuffHeader::MakePanel(wxScrolledWindow *panel)
{
	wxSizer *sizer = new wxBoxSizer(wxVERTICAL);

	if (this->prefix == "Behavior:") {
		// Play/Pause button
		wxButton *playbutton = new wxButton(panel,ID_PLAYPARTICLES,p->play ? "Pause" : "Play");
		playbutton->SetClientData(p);
		sizer->Add(playbutton,0,wxEXPAND,0);
	} else if (this->prefix == "Shader:") {
		// Show/NoShow button
		wxButton *showbutton = new wxButton(panel,ID_SHOWPARTICLES,p->show ? "Hide" : "Show");
		showbutton->SetClientData(p);
		sizer->Add(showbutton,0,wxEXPAND,0);
	}

	// Create a choicebox of names and the pointer p into which the new particlestuff should be inserted
	ParticleStuffRegistry *stuffs = &PARTICLESTUFF_REGISTRY;
	wxChoice *stuffchoice = new wxChoice(panel, ID_NEWPARTICLESTUFF);
	for (ParticleStuffRegistry::iterator stuff = stuffs->begin();stuff != stuffs->end() ; stuff++) {
		wxString stuffname = stuff->first.c_str();
		if (stuffname.StartsWith(prefix))
			stuffchoice->Append(stuff->first.data(),p);
	}

	wxSizer *rowsizer = new wxBoxSizer(wxHORIZONTAL);
	rowsizer->Add(new wxStaticText(panel,-1,_T("New ") + prefix),0,0,0);
	rowsizer->Add(stuffchoice,1,0,0);

	sizer->Add(rowsizer,0,wxEXPAND,0);

	panel->SetSizer(sizer);
	panel->Layout();
}

void ParticleWidget::OnPlayParticles(wxCommandEvent &event)
{
	wxButton *button = dynamic_cast<wxButton *>(event.GetEventObject());
	if (!button) return;
	Particles *p = (Particles *)button->GetClientData();
	p->play = !p->play;
	button->SetLabel(p->play ? "Pause" : "Play");
}

void ParticleWidget::OnShowParticles(wxCommandEvent &event)
{
	wxButton *button = dynamic_cast<wxButton *>(event.GetEventObject());
	if (!button) return;
	Particles *p = (Particles *)button->GetClientData();
	p->show = !p->show;
	button->SetLabel(p->show ? "Hide" : "Show");
}

void ParticleTreeParamButton::OnParticleTreeParamButton(wxCommandEvent &event)
{
	if (cb) cb->onbuttonpress();
}

void ParticleTreeParamComboBox::OnParticleTreeComboBoxSelectingItem(wxCommandEvent &event)
{
	std::string selection(GetValue().ToStdString());
	if (cb) cb->itemselectFromInterface(selection);
}

void
ParticleStuffItemData::MakePanel(wxScrolledWindow *panel)
{
	wxFlexGridSizer *flexsizer = new wxFlexGridSizer(0,2,1,10);
	flexsizer->AddGrowableCol(0);
	flexsizer->AddGrowableCol(1);

	for (unsigned int i = 0; i < stuff->params.size(); i++) {
		wxStaticText *label = new wxStaticText(panel,-1,stuff->params[i]->name.c_str(),wxDefaultPosition,wxDefaultSize,wxALIGN_RIGHT);
		label->SetToolTip(stuff->params[i]->desc.c_str());
		flexsizer->Add(label,1,wxEXPAND);

		// Would have preferred to do this with an overloaded funcion based on the type of stuff->params[i] but not supported by C++.
		// Have to resort to ugly dynamic_cast because PSParam should not have any reference to the widget library accessing it. -jch
		if (PSParamButton *parambutton = dynamic_cast<PSParamButton *>(stuff->params[i])) {
			ParticleTreeParamButton *ptpb = new ParticleTreeParamButton(panel,parambutton->shortname,parambutton->desc);
			ptpb->cb = parambutton->func;
			flexsizer->Add(ptpb,1,wxEXPAND);
		} else if (PSParamBool *parambool = dynamic_cast<PSParamBool *>(stuff->params[i])) {
			ParticleTreeParamCheckBox *ptpcb = new ParticleTreeParamCheckBox(panel,parambool);
			ptpcb->SetToolTip(parambool->desc);
			flexsizer->Add(ptpcb,1,wxEXPAND);
		} 
		else if (PSParamComboBox * paramCombo = dynamic_cast<PSParamComboBox *>(stuff->params[i]))
		{				
			wxArrayString choices;
			for (unsigned int k=0;k<paramCombo->func->nbChoices();++k)
			{
				choices.Add(paramCombo->func->getChoiceString(k));
			}
		
			ParticleTreeParamComboBox *ptcb = new ParticleTreeParamComboBox(panel, choices,
																			paramCombo->shortname, paramCombo->desc);
			ptcb->cb = paramCombo->func;
			//we set the interface to the selected value
			//this is important if we leave and come back, or if the state is loaded
			ptcb->SetValue(paramCombo->func->getSelectedString());
			flexsizer->Add(ptcb,1,wxEXPAND);
		}
		else {
			ParticleTreeStuffParam *ptsp = new ParticleTreeStuffParam(panel,stuff->params[i]);
			ptsp->SetToolTip(stuff->params[i]->desc);
			flexsizer->Add(ptsp,1,wxEXPAND);
		}
	}

	for (unsigned int i = 0; i < stuff->attachedattributes.size(); i++) {
		wxStaticText *label = new wxStaticText(panel,-1,stuff->attachedattributes[i]->name.c_str(),wxDefaultPosition,wxDefaultSize,wxALIGN_RIGHT);
		label->SetToolTip(stuff->attachedattributes[i]->desc.c_str());
		flexsizer->Add(label,1,wxEXPAND);
		ParticleTreeStuffAA *ptsaa = new ParticleTreeStuffAA(panel,stuff->attachedattributes[i]);
		ptsaa->SetToolTip(stuff->attachedattributes[i]->desc.c_str());
		flexsizer->Add(ptsaa,1,wxEXPAND);
	}

	int selected = stuff->ps->selectedParticle;

	if (stuff->perparticle.size()) {
		flexsizer->Add(new wxStaticText(panel,-1,_T("Selected Particle:"),wxDefaultPosition,wxDefaultSize,wxALIGN_RIGHT),1,wxEXPAND);
		flexsizer->Add(new ParticleTreeShowSelected(panel,stuff->ps),1,wxEXPAND);
		for (unsigned int i = 0; i < stuff->perparticle.size(); i++) {
			wxStaticText *label = new wxStaticText(panel,-1,stuff->perparticle[i]->name.c_str(),wxDefaultPosition,wxDefaultSize,wxALIGN_RIGHT);
			label->SetToolTip(stuff->perparticle[i]->desc.c_str());
			flexsizer->Add(label,1,wxEXPAND);
			ParticleTreeStuffParamPerParticle *ptsppp = new ParticleTreeStuffParamPerParticle(panel,stuff->perparticle[i],selected);
			ptsppp->SetToolTip(stuff->perparticle[i]->desc.c_str());
			flexsizer->Add(ptsppp,1,wxEXPAND);
		}
	}

	panel->SetSizer(flexsizer);
	panel->Layout();
}

void
ImplicitInterrogatorItemData::MakePanel(wxScrolledWindow *panel)
{
	DoubleVector q;
	NameVector qnames;
	wxString buf;

	wxSizer *sizer = new wxBoxSizer(wxVERTICAL);
	wxSizer *rowsizer;

	stuff->getq(q);

	Surfaces *surfaces = &(wxGetApp().surfaces);
	wxChoice *impchoice = new wxChoice(panel,ID_IMPINTCHANGED);
	for (unsigned int i = 0; i < (*surfaces).size(); i++) {
		impchoice->Append((*surfaces)[i]->getObjectName().c_str(), (void *)stuff);
	}
	impchoice->SetSelection((int)q[0]);

	rowsizer = new wxBoxSizer(wxHORIZONTAL);
	rowsizer->Add(new wxStaticText(panel,-1,_T("Implicit:")),1,wxEXPAND,0);

	rowsizer->Add(impchoice);
	sizer->Add(rowsizer,0,wxEXPAND,0);

	panel->SetSizer(sizer);
	panel->Layout();
}

void
ParticleWidget::OnImpIntChanged(wxCommandEvent &event)
{
	ImplicitInterrogator *ii = (ImplicitInterrogator *)event.GetClientData();
	double impindex;
	Surfaces *surfaces = &(wxGetApp().surfaces);

	ii->getq(&impindex);
	impindex = event.GetSelection();
	ii->setq(&impindex);
	ii->setImplicit((Implicit *)((*surfaces)[(int) impindex]));
}

void
ParticleWidget::OnNewParticleSystem(wxCommandEvent &event)
{
	ParticleSystem *ps = new ParticleSystem();
	ps->name = event.GetString().c_str();
	ps->surfaces = &(wxGetApp().surfaces);
	ParticleSystems *pss = &(wxGetApp().particlesystems);
	pss->push_back(ps);
	pss->surfacesAttach();
	wxTreeItemId _id = partree->GetRootItem();
	partree->Add(_id,ps);
	partree->Refresh();
//	Reload();
}

void
ParticleWidget::OnNewParticles(wxCommandEvent &event)
{
	ParticleSystemItemData *psit = (ParticleSystemItemData *)(partree->GetItemData(partree->GetSelection()));
	if (!psit) return;

	wxString name = event.GetString();

	Particles *p = new Particles(psit->ps,std::string(name.GetData()));

	wxTreeItemId _id = partree->GetSelection();
	partree->Add(_id,p);
	partree->Refresh();
//	Reload();
}

void
ParticleWidget::OnCopyParticles(wxCommandEvent &event)
{
	ParticlesItemData *pit = (ParticlesItemData *)(partree->GetItemData(partree->GetSelection()));
	if (!pit || !pit->p) return;

	std::string src_name = event.GetString().ToStdString();
	Particles *p = pit->p->particleSystem->findParticles(src_name);
	pit->p->copyFrom(p);

	// Add newly copied attributes, behaviors and shaders to existing particle tree

	// Currently selected node must be the particles node
	ParticleTreeItem pnode = partree->GetSelection();

	// It must have only three children which are the title nodes for
	// attributes, behaviors and shaders.
	wxTreeItemIdValue cookie;
	ParticleTreeItem attrnode = partree->GetFirstChild(pnode,cookie);
	ParticleTreeItem behanode = partree->GetNextChild(pnode,cookie);
	ParticleTreeItem shadnode = partree->GetNextChild(pnode,cookie);

	partree->Add(attrnode,pit->p->attributes);
	partree->Add(behanode,pit->p->behaviors);
	partree->Add(shadnode,pit->p->shaders);

	// Show the particles node with its new attr, behaviors and shaders below it
	partree->EnsureVisible(pnode);
}

void
ParticleWidget::OnNewParticleStuff(wxCommandEvent &event)
{
	Particles *p = (Particles *)event.GetClientData();

	wxString stuffname = event.GetString();
	std::auto_ptr<ParticleStuff> temp = NEW_PARTICLESTUFF(stuffname.ToStdString());
	ParticleStuff *stuff = temp.release();
	stuff->setParticleSystem(p);

	/* Warning: this only works when the current selected tree node
	 * is an "Attributes", "Behaviors" or "Shaders" header node.
	 */
	wxTreeItemId _id = partree->GetSelection(); //added for sutpid gcc --ms
	partree->Add(_id,stuff);

	p->attachAttributes();

	// new attributes could have been created
	// this code clears and refills the attributes portion of the particle tree
	ParticleTreeItem pnode = partree->GetItemParent(partree->GetSelection());
	wxTreeItemIdValue cookie;
	ParticleTreeItem attrnode = partree->GetFirstChild(pnode,cookie);
	partree->DeleteChildren(attrnode);
	partree->Add(attrnode,p->attributes);

	partree->Expand(partree->GetSelection());
}
