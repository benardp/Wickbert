#include "wxModeler.h"
#include "LogFrame.h"
#include "SurfaceTree.h"
#include "Surface/Implicit/Implicit.h"
#include "Surface/OpenMesh/SurfaceMesh.h"

BEGIN_EVENT_TABLE(SurfaceTree, wxSplitterWindow)
	EVT_LIST_END_LABEL_EDIT(ID_SURFACETREE,SurfaceTree::OnRename)
	EVT_LIST_ITEM_SELECTED(ID_SURFACETREE,SurfaceTree::OnSelect)
	EVT_CHOICE(ID_NEWSURFACE,SurfaceTree::OnNew)
	EVT_LIST_KEY_DOWN(ID_SURFACETREE,SurfaceTree::OnKeyPress)
	EVT_BUTTON(ID_UPDEGREE,SurfaceTree::OnUpDegree)
	EVT_BUTTON(ID_DOWNDEGREE,SurfaceTree::OnDownDegree)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(SurfaceTreeParamButton, wxButton)
EVT_BUTTON(ID_SURFACETREEPARAMBUTTON,SurfaceTreeParamButton::OnSurfaceTreeParamButton)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(SurfaceTreeParamComboBox, wxComboBox)
EVT_COMBOBOX(ID_SURFACETREEPARAMCOMBOBOX, SurfaceTreeParamComboBox::OnSurfaceTreeComboBoxSelectingItem)
END_EVENT_TABLE()



void
SurfaceTree::OnKeyPress(wxListEvent &event)
{
		long id = event.GetIndex();
		Surface *s = (Surface *)(surflist->GetItemData(id));
		Surfaces *surfaces = &(wxGetApp().surfaces);
		/* This should be a member function of Surfaces */
		Surfaces::iterator sit = find(surfaces->begin(),surfaces->end(),s);
	
	switch (event.GetKeyCode()) {

		case WXK_DELETE:
		//added for wxWidgets 2.6.3 (I don't know why 
		//Delete key generates event with key code
		// 389 instead of WXK_DELETE) ---ms
		case 389:
			if (sit == surfaces->end()) {
				wxMessageBox("Currently selected surface not found.\n"
							 "Removing from displayed list anyway.","Delete Error",wxOK | wxICON_HAND);
			} else {
				surfaces->erase(sit);
			}

			/* erase even if not found because then it shouldn't be there anyway. */
			surflist->DeleteItem(id);
			surfaces->attach();
			/* should probably reattach all the particlesystem stuff too. */
		break;
	}
}

void
SurfaceTree::OnRename(wxListEvent &event)
{
	long id = event.GetIndex();
	Surface *s = (Surface *)(surflist->GetItemData(id));
	wxString newname = event.GetLabel();
	s->setObjectName(newname.ToStdString());
}

void
SurfaceTree::OnNew(wxCommandEvent &event)
{
	wxString surfname = event.GetString();

	Surfaces *surfaces = &(wxGetApp().surfaces);

	std::auto_ptr<Implicit> imp = NEW_IMPLICIT(surfname.ToStdString());
	surfaces->push_back(imp.release());

	surflist->DeleteAllItems();
	Add(surfaces);
}

SurfaceTree::SurfaceTree(wxWindow *parent, wxWindowID id, wxPoint pos, wxSize size)
: wxSplitterWindow(parent,-1,wxDefaultPosition,wxDefaultSize,wxSP_3D)
{
	ImplicitRegistry *impReg = &IMPLICIT_REGISTRY;

	selected = NULL;

	SetBackgroundColour(*wxWHITE);
	surflistpan = new wxScrolledWindow(this);
	surfparpan =  new wxScrolledWindow(this);

	wxBoxSizer *slpsizer = new wxBoxSizer(wxVERTICAL);
	
/* Don't need this now that there is a Surface->New menu */
#ifdef NEWSURFACE_CHOICEBOX
	wxBoxSizer *rowsizer = new wxBoxSizer(wxHORIZONTAL);
	rowsizer->Add(new wxStaticText(surflistpan,-1,"New... "),0,0,0);
	newsurf = new wxChoice(surflistpan, ID_NEWSURFACE, wxDefaultPosition, wxSize(200,-1) );
	for (ImplicitRegistry::iterator impobj = impReg->begin();impobj != impReg->end() ; impobj++) {
		std::auto_ptr<Implicit> temp = NEW_IMPLICIT(impobj->first);
		newsurf->Append(impobj->first.data(),temp.release());
	}
	rowsizer->Add(newsurf,1,0,0);
	slpsizer->Add(rowsizer,0,wxEXPAND,0);
#endif

	surflist = new wxListCtrl(surflistpan,ID_SURFACETREE,wxDefaultPosition,wxDefaultSize,wxLC_REPORT | wxLC_EDIT_LABELS);
	surflist->InsertColumn(0, _T("Name"));
    surflist->InsertColumn(1, _T("Type"));
	slpsizer->Add(surflist,1,wxEXPAND,0);

	Surfaces *surfaces = &(wxGetApp().surfaces);
	Add(surfaces);

	surflistpan->SetSizerAndFit(slpsizer);
	
	SetMinimumPaneSize(20);
	SplitHorizontally(surflistpan,surfparpan);
}

void
SurfaceTree::Add(Surfaces *surfs)
{
	for (unsigned int i = 0; i < surfs->size(); i++) {
//		surfitem = new wxListItem();
//		surfitem->SetData(((*surfs)[i]));
//		surfitem->SetText((*surfs)[i]->name().c_str());
		Add(surfs->at(i));
	}
}

void
SurfaceTree::Add(Surface *s)
{
	long id = surflist->InsertItem(surflist->GetItemCount(),s->getObjectName().c_str());
	surflist->SetItemData(id,(long)s);

	wxListItem item;
	item.m_mask = wxLIST_MASK_TEXT;
	item.m_text = s->name().c_str();
	item.m_col = 1;
	item.m_itemId = id;

	surflist->SetItem(item);
}

void
SurfaceTree::OnSelect(wxListEvent &event)
{
	surfparpan->DestroyChildren();

	long id = event.GetIndex();
	selected = (Surface *)(surflist->GetItemData(id));

	DrawParamPanel(selected);

	Refresh();
}

void
SurfaceTree::DrawParamPanel(Surface *s)
{
	Implicit *imp = dynamic_cast<Implicit *>(s);
	SurfaceMesh *mesh = dynamic_cast<SurfaceMesh *>(s);

	wxSizer *sizer = new wxBoxSizer(wxVERTICAL);
	wxSizer *rowsizer;

	/* Compute max label width */
	int width = 0;
	wxStaticText *test = new wxStaticText(surfparpan,-1,"",wxDefaultPosition);

	for (unsigned int i = 0; i < s->params.size(); i++) {
		int w,h;
		test->SetLabel(s->params[i]->name.c_str());
		test->GetSize(&w,&h);
		if (w > width) width = w;
	}
	if (imp ) {
		NameVector qnames;
		imp->getqname(qnames);
		
		for (unsigned int i = 0; i < qnames.size(); i++) {
			int w,h;
			test->SetLabel(qnames[i].c_str());
			test->GetSize(&w,&h);
			if (w > width) width = w;
		}
	}
	test->Destroy();

	/* Now draw the panel */
	for (unsigned int i = 0; i < s->params.size(); i++) {		
		rowsizer = new wxBoxSizer(wxHORIZONTAL);
		rowsizer->Add(new wxStaticText(surfparpan,-1,s->params[i]->name.c_str(),wxDefaultPosition,wxSize(width,-1),wxALIGN_RIGHT),0,0,0);
		if (SurfParamButton *parambutton = dynamic_cast<SurfParamButton *>(s->params[i])) {
			SurfaceTreeParamButton *ptpb = new SurfaceTreeParamButton(surfparpan,parambutton->shortname,parambutton->desc);
			ptpb->cb = parambutton->func;
			rowsizer->Add(ptpb);
		}
		else if (SurfParamComboBox * paramCombo = dynamic_cast<SurfParamComboBox*>(s->params[i]))
		{				
			wxArrayString choices;
			for (unsigned int i=0;i<paramCombo->func->nbChoices();++i)
			{
				choices.Add(paramCombo->func->getChoiceString(i));
			}
		
			SurfaceTreeParamComboBox *ptcb = new SurfaceTreeParamComboBox(	surfparpan, choices,
																			paramCombo->shortname, paramCombo->desc);
			ptcb->cb = paramCombo->func;
			//we set the interface to the selected value
			//this is important if we leave and come back, or if the state is loaded
			ptcb->SetValue(paramCombo->func->getSelectedString());
			rowsizer->Add(ptcb);
		}
		else
		{
			rowsizer->Add(new SurfParamTextCtrl(surfparpan,wxString(s->params[i]->get().c_str()),s,i));
		}
		sizer->Add(rowsizer,0,wxEXPAND,0);
	}

	if (imp) {

		Algebraic *alg = dynamic_cast<Algebraic *>(imp);
		if (alg) {
			rowsizer = new wxBoxSizer(wxHORIZONTAL);
			wxString deglabel;
			deglabel << "Degree: " << alg->degree() << " ";
			rowsizer->Add(new wxStaticText(surfparpan,-1,deglabel),0,0,0);
			wxButton *upbut = new wxButton(surfparpan,ID_UPDEGREE,"+",wxDefaultPosition,wxSize(30,30));
			upbut->SetClientData(alg);
			rowsizer->Add(upbut,0,0,0);
			wxButton *downbut = new wxButton(surfparpan,ID_DOWNDEGREE,"-",wxDefaultPosition,wxSize(30,30));
			downbut->SetClientData(alg);
			rowsizer->Add(downbut,0,0,0);
			sizer->Add(rowsizer,0,wxEXPAND,0);
		}

		DoubleVector q;
		NameVector qnames;
		wxString buf;

		imp->getq(q);
		imp->getqname(qnames);

		for (unsigned int i = 0; i < imp->qlen(); i++) {
			rowsizer = new wxBoxSizer(wxHORIZONTAL);
			rowsizer->Add(new wxStaticText(surfparpan,-1,qnames[i].c_str(),wxDefaultPosition,wxSize(width,-1),wxALIGN_RIGHT),0,0,0);
			buf.Printf("%g",q[i]);
			rowsizer->Add(10,10);
			rowsizer->Add(new ImpParamTextCtrl(surfparpan,buf,imp,i));
			sizer->Add(rowsizer,0,wxEXPAND,0);
		}
 
	} else if (mesh) {
		wxString verts; verts << "Vertices: " << mesh->n_vertices();
		sizer->Add(new wxStaticText(surfparpan,-1,verts),0,wxLEFT | wxTOP,10);
		wxString faces; faces << "Faces: " << mesh->n_faces();
		sizer->Add(new wxStaticText(surfparpan,-1,faces),0,wxLEFT | wxTOP,10);
		wxString edges; edges << "Edges: " << mesh->n_edges();
		sizer->Add(new wxStaticText(surfparpan,-1,edges),0,wxLEFT | wxTOP,10);
	}
	surfparpan->SetSizer(sizer);
	surfparpan->SetScrollbars(1,1,width+210,1000);
	surfparpan->Layout();
}

BEGIN_EVENT_TABLE(ImpParamTextCtrl, wxTextCtrl)
	EVT_TEXT_ENTER(ID_IMPLICITPARAM,ImpParamTextCtrl::OnValue)
// need to fix for 2.6.0
	EVT_UPDATE_UI(wxID_ANY,ImpParamTextCtrl::OnUpdateUI)
END_EVENT_TABLE()

void
ImpParamTextCtrl::OnValue(wxCommandEvent &event)
{
	ImpParamTextCtrl *sptc = (ImpParamTextCtrl *)(event.GetEventObject());

	DoubleVector q;

	sptc->imp->getq(q);
	sptc->GetValue().ToDouble(&(q[sptc->i]));
	sptc->imp->setq(q);
}

void
ImpParamTextCtrl::OnUpdateUI(wxUpdateUIEvent &event)
{
	if (FindFocus() == this) return;

	ImpParamTextCtrl *sptc = (ImpParamTextCtrl *)(event.GetEventObject());

	DoubleVector q;

	sptc->imp->getq(q);

	// Dunno why this is necessary, but I get an error otherwise
	if ((unsigned int)sptc->i >= q.size())
		return;

	wxString buf;
	sptc->SetValue(buf.Format("%g",q[sptc->i]));
}

BEGIN_EVENT_TABLE(SurfParamTextCtrl, wxTextCtrl)
	EVT_TEXT_ENTER(ID_SURFACEPARAM,SurfParamTextCtrl::OnValue)
END_EVENT_TABLE()

void
SurfParamTextCtrl::OnValue(wxCommandEvent &event)
{
	SurfParamTextCtrl *sptc = (SurfParamTextCtrl *)(event.GetEventObject());

	sptc->s->params[sptc->i]->set(sptc->GetValue().ToStdString());

	Surfaces *surfaces = &(wxGetApp().surfaces);
	surfaces->attach();
	//now we wanna see whether the attachments were successful, thus we need to update the value 
	//in the text output.
	SurfStringParam * surfString=dynamic_cast<SurfStringParam*> (sptc->s->params[sptc->i]);
	if (surfString!=0)
	{
		sptc->SetValue(surfString->v);
	}
}

void SurfaceTreeParamButton::OnSurfaceTreeParamButton(wxCommandEvent &event)
{
	if (cb) cb->onbuttonpress();
}

void SurfaceTreeParamComboBox::OnSurfaceTreeComboBoxSelectingItem(wxCommandEvent &event)
{
	std::string selection(GetValue().ToStdString());
	if (cb) cb->itemselectFromInterface(selection);
}
	