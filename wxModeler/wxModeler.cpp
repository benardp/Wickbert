#include "wxModeler.h"
#include "LogFrame.h"

#ifdef WB_USE_CGAL
#include "ClipartFrame.h"
#endif

#include "DisplayWindow.h"
#include "ParticleTree.h"
#include "SurfaceTree.h"
#include "Surface/Implicit/ImpFileManager.h"

#ifndef _WIN32
#undef check
#endif

#include "Surface/OpenMesh/SurfaceMesh.h"

#include <iostream>

IMPLEMENT_APP(Modeler)

BEGIN_EVENT_TABLE(Modeler, wxApp)

END_EVENT_TABLE()

BEGIN_EVENT_TABLE(MainFrame, wxFrame)
	EVT_MENU(wxID_OPEN, MainFrame::OnOpen)
	EVT_MENU(wxID_SAVE, MainFrame::OnSave)
	EVT_MENU(ID_SAVE_SVG, MainFrame::OnSaveSVG)
	EVT_MENU(ID_SAVE_GL, MainFrame::OnSaveGL)
	EVT_MENU(ID_START_SAVE_SVG_ANIMATION, MainFrame::OnStartSaveSVGAnimation)
	EVT_MENU(ID_STOP_SAVE_SVG_ANIMATION, MainFrame::OnStopSaveSVGAnimation)
	EVT_MENU(ID_START_SAVE_GL_ANIMATION, MainFrame::OnStartSaveGLAnimation)
	EVT_MENU(ID_STOP_SAVE_GL_ANIMATION, MainFrame::OnStopSaveGLAnimation)
	EVT_MENU_RANGE(ID_NEWSURFACEMENU,ID_NEWSURFACEMENU_LAST,MainFrame::OnNewSurfaceMenu)
	EVT_MENU(ID_PAUSEALL, MainFrame::OnPauseAll)
	EVT_MENU(ID_PLAYALL, MainFrame::OnPlayAll)
	EVT_MENU(ID_REFRESHPARWIN, MainFrame::OnRefreshParWin)
	EVT_MENU(ID_PARPROFILE, MainFrame::OnParProfile)
	EVT_MENU(ID_OPENGLEXTS, MainFrame::OnOpenGLExts)
	EVT_MENU(ID_ZOOMIN, MainFrame::OnZoomIn)
	EVT_MENU(ID_ZOOMOUT, MainFrame::OnZoomOut)
	EVT_MENU(wxID_EXIT, MainFrame::OnQuit)
	EVT_MENU(wxID_ABOUT, MainFrame::OnAbout)
    EVT_IDLE(MainFrame::OnIdle)
	EVT_MENU(ID_SHOWLOG, MainFrame::OnShowLog)
	EVT_MENU(ID_SHOWCLIPART, MainFrame::OnShowClipArt)
END_EVENT_TABLE()

////////////////////////////////////////
//Modeler implementation
////////////////////////////////////////

bool Modeler::OnInit()
{
	/* initialize the factories */
	impReg = &IMPLICIT_REGISTRY;
	psReg = &PARTICLESTUFF_REGISTRY;

	surfaces.setParticleSystems(&particlesystems);

	mainframe = new MainFrame("wxModeler",wxDefaultPosition,wxSize(1000,600));
	mainframe->Show(TRUE);
    SetTopWindow(mainframe);

	if (argc == 2) {
		wxString filename = argv[1];
		if (filename.Right(2) == "wb")
			mainframe->LoadWB(filename);
	} else {
		/* build a default implicit */
		std::auto_ptr<Implicit> temp=NEW_IMPLICIT("Algebraic:Torus");
		Implicit *imp=temp.release();
		surfaces.push_back(imp);
		mainframe->impwin->Add(&surfaces);

		mainframe->LoadPar("defaultparticlesystem.par");
	}
	logWindow = NULL;

#ifdef WB_USE_CGAL
	clipArtFrame = NULL;
#endif

	captureSVGflag = false;
	captureCanvasFlag = false;
	return TRUE;
}

int Modeler::OnExit()
{
	//here delete objects in the future
	return 0;
}

////////////////////////////////////////////////
//MainFrame Implementation
////////////////////////////////////////////////

/** Constructs a cascading menu of new surface objects.
 */

wxMenu *
MainFrame::NewSurfaceMenu()
{
	wxMenu *newsurfmenu = new wxMenu;
	ImplicitRegistry *impReg = &IMPLICIT_REGISTRY;

	//idtosurfname.reserve(ID_NEWSURFACEMENU_LAST - ID_NEWSURFACEMENU); 
	idtosurfname.assign(ID_NEWSURFACEMENU_LAST - ID_NEWSURFACEMENU, NULL); //--reserve is not enough for line 146--ms

	int uniqueid = ID_NEWSURFACEMENU + 1;
	for (ImplicitRegistry::iterator impobj = impReg->begin();impobj != impReg->end() ; impobj++) {
		std::auto_ptr<Implicit> temp = NEW_IMPLICIT(impobj->first);
		temp.release();
		std::string name = impobj->first;
		std::string pre;
		wxMenu *currentmenu = newsurfmenu;
		bool done = false;
		while (!done) {
			/* separate name into pre:name */
			int colon = name.find(':');
			if (colon == -1) {
				done = true;
				pre = name;
				name = "";
			} else {
				pre = name.substr(0,colon+1);
				name = name.substr(colon+1);
			}

			int id = currentmenu->FindItem(pre.c_str());
			if (id == wxNOT_FOUND) {
				id = uniqueid++;
				if (done) {
					currentmenu->Append(id,pre.c_str());
					idtosurfname[id - ID_NEWSURFACEMENU] = new std::string(impobj->first);
				} else {
					wxMenu *submenu = new wxMenu;
					currentmenu->Append(id,pre.c_str(),submenu);
				}
			}
			if (!done) {
				wxMenuItem *item = currentmenu->FindItem(id);
				currentmenu = item->GetSubMenu();
			}
		}			
	}

	return newsurfmenu;
}

void
MainFrame::OnNewSurfaceMenu(wxCommandEvent &event)
{
	int id = event.GetId();
	std::string *surfname = idtosurfname[id - ID_NEWSURFACEMENU];

	Surfaces *surfaces = &(wxGetApp().surfaces);

	std::auto_ptr<Implicit> imp = NEW_IMPLICIT(surfname->c_str());
	surfaces->push_back(imp.release());

#if 0
	impwin->surflist->DeleteAllItems();
	impwin->Add(surfaces);
#else
	impwin->Add(surfaces->at(surfaces->size()-1));
#endif
}

MainFrame::MainFrame(const wxString &title, const wxPoint &pos, const wxSize &size)
: wxFrame(NULL, -1, title, pos, size)
{
	wxMenuBar *menubar = new wxMenuBar;

	wxMenu *filemenu = new wxMenu;
	filemenu->Append(wxID_OPEN, "&Open...");
	filemenu->Append(wxID_SAVE, "&Save...");
	filemenu->Append(ID_SAVE_GL, "Save GL Canvas");
	filemenu->Append(ID_SAVE_SVG, "Save SVG");
	filemenu->Append(ID_START_SAVE_GL_ANIMATION, "Start Save GL Animation");
	filemenu->Append(ID_STOP_SAVE_GL_ANIMATION, "Stop Save GL Animation");
	filemenu->Append(ID_START_SAVE_SVG_ANIMATION, "Start Save SVG Animation");
	filemenu->Append(ID_STOP_SAVE_SVG_ANIMATION, "Stop Save SVG Animation");
	filemenu->Append(wxID_ABOUT, "&About...");
	filemenu->Append(wxID_EXIT, "E&xit");

	menubar->Append(filemenu,"&File");

	wxMenu *viewmenu = new wxMenu;
	viewmenu->Append(ID_ZOOMIN, "Zoom &in");
	viewmenu->Append(ID_ZOOMOUT, "Zoom &out");
	viewmenu->Append(ID_OPENGLEXTS,"OpenGL &Extensions");
	viewmenu->Append(ID_SHOWLOG,"Show &Log Window");
	viewmenu->Append(ID_SHOWCLIPART,"Show &ClipArt Window");
	
	menubar->Append(viewmenu,"&View");

	wxMenu *parmenu = new wxMenu;
	parmenu->Append(ID_PAUSEALL, "Pause All");
	parmenu->Append(ID_PLAYALL, "Play All");
	parmenu->Append(ID_REFRESHPARWIN, "Refresh Tree");
	parmenu->Append(ID_PARPROFILE,"&Profile...");

	menubar->Append(parmenu,"&Particles");

	wxMenu *surfmenu = new wxMenu;
	surfmenu->Append(ID_NEWSURFACEMENU,"New",NewSurfaceMenu());

	menubar->Append(surfmenu,"&Surfaces");

	SetMenuBar(menubar);

	/* three-way vertical split */
	wxSplitterWindow *glparimp = new wxSplitterWindow(this,-1,wxDefaultPosition,wxDefaultSize,wxSP_3D);
	displaywin = new DisplayWindow(glparimp,-1,wxDefaultPosition,wxDefaultSize);
	wxSplitterWindow *parimp = new wxSplitterWindow(glparimp,-1,wxDefaultPosition,wxDefaultSize,wxSP_3D);

	glparimp->SplitVertically(displaywin,parimp);
	parwin = new ParticleWidget(parimp);
	impwin = new SurfaceTree(parimp);
	parimp->SplitVertically(parwin,impwin);
	
	//spliting order seems to matter for wxMac, so split horizontally after the vertical split
	parwin->SplitHorizontally(parwin->partree,parwin->parpanel); //this was in the parwin constructor
		
#if 0
	wxToolBar *toolbar = CreateToolBar(wxHORIZONTAL, ID_TOOLBAR);
	wxComboBox *combo = new wxComboBox(toolbar, ID_COMBO, _T(""), wxDefaultPosition, wxSize(200,-1) );
	for (ImplicitRegistry::iterator impobj = impReg->begin();impobj != impReg->end() ; impobj++) {
		std::auto_ptr<Implicit> temp = NEW_IMPLICIT(impobj->first);
		combo->Append(impobj->first.data(),temp.release());
	}
	toolbar->AddControl(combo);
#endif

	displaywin->SetupView(FALSE,0,0);
	
	CreateStatusBar();
	SetStatusText("Everything's going to be just fine.");
}

void MainFrame::OnOpen(wxCommandEvent &event)
{
	ParticleSystems *particlesystems = &(wxGetApp().particlesystems);
	Surfaces *surfaces = &(wxGetApp().surfaces);
	//WHY??? I comment this out for the moment. This might lead to a crash, when no implicit is present and we load something... - Elmar 
	//Implicit *imp = dynamic_cast<Implicit *>((*surfaces)[0]);

	wxFileDialog *filedialog = new wxFileDialog(this, "Open a File", "", "",
		"Particles + Implicits (*.wb)|*.wb|"
		"Particle Systems (*.par)|*.par|"
		"Implicit Surfaces (*.imp)|*.imp|"
		"Meshes (*.obj)|*.obj|"
		"Other files (*.*)|*.*", wxFD_OPEN );

	if (filedialog->ShowModal() == wxID_CANCEL)
		return;

	wxString filename = filedialog->GetPath();

	/* Redirect cerr into the string errors */
	std::stringbuf errors;
	std::streambuf *cerrbuforig = std::cerr.rdbuf();
	std::cerr.rdbuf(&errors);

	if (filename.Right(2) == "wb") {
		LoadWB(filename);
	} else if (filename.Right(3) == "par") {
		LoadPar(filename);
	} else if (filename.Right(3) == "imp") {
		/* load an implicit */
		Surfaces newsurfs;
        ImpFileManager reader;
		bool check=reader.readFile(filename.ToStdString(), newsurfs);
		for (Surfaces::iterator surf_it = newsurfs.begin(); surf_it != newsurfs.end(); surf_it++) {
			surfaces->push_back(*surf_it);
		}
		impwin->Add(&newsurfs);
#if 0
        Surfaces newsurfs;
		std::ifstream in(filename);
		while (in.peek() != EOF)
			in >> &newsurfs;
		impwin->Add(&newsurfs);
		surfaces->insert(surfaces->end(),newsurfs[0]);
#endif
	} else if (filename.Right(3) == "obj") {
		/* load a meshed surface */
		Surface *s = new SurfaceMesh(filename.ToStdString());
		surfaces->push_back(s);
		impwin->Add(s);
	}

	std::cerr.flush();

	if (errors.str().size()) {
		wxMessageBox((errors.str()).c_str(),"Error Loading File",wxICON_ERROR | wxOK,this);
	}

	/* Restore cerr */
	std::cerr.rdbuf(cerrbuforig);
}

void
MainFrame::LoadWB(wxString filename)
{
	ParticleSystems *particlesystems = &(wxGetApp().particlesystems);
	Surfaces *surfaces = &(wxGetApp().surfaces);

	std::ifstream in(filename.ToStdString());
	std::string token;
	if (!(in >> token) || token != "Wickbert-0.3") {
		wxMessageBox("Wrong version.\n"
			"Expected \"Wickbert-0.3\" but read " + token + ".");
	} else {
		while (in >> token) {
			if (token == "particlesystem") {
				ParticleSystem *newps = new ParticleSystem(particlesystems);
				/* need to read in name of particle system */
				char c;
				if (in >> c && c == '"') {
					in >> std::noskipws;
					token = "";
					while (in >> c && c != '"')
						token += c;
					newps->name = token;
					in >> std::skipws;
				}
				if (in >> c && c != '{') {
					std::cerr << "Expected '{' after particle system name";
					break;
				}
				if (in >> newps) {
					particlesystems->push_back(newps);
					//ParticleTreeItem root = parwin->partree->GetRootItem();
					//parwin->partree->Add(root,newps);
				}
			} else if (token == "surface") {
				Surfaces newsurfs;
				in >> &newsurfs;		/* reads just one */
				if (newsurfs.size()) {
					surfaces->insert(surfaces->end(),newsurfs[0]);
					impwin->Add(&newsurfs);
				}
			}else if( token == "mesh" ){
				// add code here to load non-implicit surfaces
				SurfaceMesh *mesh = new SurfaceMesh;
				if ( ! mesh->readStream(in) )
				{
					// couldn't find mesh?
					wxFileDialog *filedialog = new wxFileDialog(this, "Open a Mesh File", "", "",
						"Meshes (*.obj)|*.obj|"
						"Other files (*.*)|*.*", wxFD_OPEN);

					if (filedialog->ShowModal() != wxID_CANCEL) {
						wxString filename = filedialog->GetPath();
						mesh->readFile(filename.ToStdString().c_str());
					}
				}

				surfaces->insert(surfaces->end(), mesh);
				impwin->Add(mesh);
			} else {
				wxMessageBox("Read \"" + token + "\" when expecting surface or particlesystem.");
			}
		}
		//parwin->Reload();
		parwin->partree->Reload();
		particlesystems->setSurfaces(surfaces);
		particlesystems->attachAttributes();
		surfaces->attach();
	}
}

void MainFrame::LoadPar(wxString filename)
{
 	ParticleSystems *particlesystems = &(wxGetApp().particlesystems);
	ParticleSystems loaded;
	Surfaces *surfaces = &(wxGetApp().surfaces);
	std::string version;

	std::ifstream in(filename.ToStdString());
	if (!in) return;

	in >> version;
	if (version == "v0.3") {
		in >> &loaded;
	} else {
		std::cerr << "Wrong version. Expected \"v0.3\" but got \"" + version + "\" instead.";
	}
	in.close();
	ParticleSystems::iterator iter;
	if (loaded.size()) {
		//parwin->partree->Reload();
		for(iter = loaded.begin(); iter != loaded.end(); ++iter)
			(*iter)->particleSystems = particlesystems;
		parwin->partree->Add(&loaded);
		loaded.setSurfaces(surfaces);
		particlesystems->insert(particlesystems->end(),loaded.begin(),loaded.end());
		particlesystems->attachAttributes();
		
	}
}

void MainFrame::OnSave(wxCommandEvent &event)
{
	ParticleSystems *particlesystems = &(wxGetApp().particlesystems);
	Surfaces *surfaces = &(wxGetApp().surfaces);
	Implicit *imp = dynamic_cast<Implicit *>((*surfaces)[0]);

	wxFileDialog *filedialog = new wxFileDialog(this, "Save to a File", "", "",
		"Particles + Implicits (*.wb)|*.wb|"
		"Particle Systems (*.par)|*.par|"
		"Implicit Surfaces (*.imp)|*.imp|"
		"Meshed Surfaces (*.obj)|*.obj|"
		"Other files (*.*)|*.*", wxFD_SAVE);

	if (filedialog->ShowModal() == wxID_CANCEL)
		return;

	wxString filename = filedialog->GetPath();

	if (filename.Right(2) == "wb") {
		std::ofstream out(filename.ToStdString());
		if (out) {
			out << "Wickbert-0.3" << std::endl;
			for (unsigned int i = 0; i < particlesystems->size(); i++) {
				out << "particlesystem " << particlesystems->at(i);
			}
			for (unsigned int i = 0; i < surfaces->size(); i++) {
				SurfaceMesh *sm = dynamic_cast<SurfaceMesh *>(surfaces->at(i));
				if ( sm ) {
					out << "mesh ";
					sm->writeStream(out);
				}
				else
					out << "surface " << surfaces->at(i); // run some logic here to determine if this is a surface or an implicit
	//			out << "implicit " << surfaces->at(i);
			}
			out.close();
		}
	} else if (filename.Right(3) == "par") {
		std::ofstream out(filename.ToStdString());
		if (out) {
			out << "v0.3" << std::endl;
			/* save particle system */
			out << *particlesystems;
			out.close();
		}
	} else if (filename.Right(3) == "imp") {
		std::ofstream out(filename.ToStdString());
		if (out) {
			/* save surfaces */
			out << surfaces;
			out.close();
		}
	} else if (filename.Right(3) == "obj") {
		if (!impwin->selected) {
			wxMessageBox("No surface mesh selected", "Oops", wxICON_ERROR);
			return;
		}
		SurfaceMesh* mesh = dynamic_cast<SurfaceMesh *>(impwin->selected);
		if (!mesh) {
			wxMessageBox("Selected surface " + impwin->selected->name() + " not a mesh.", "Oops", wxICON_ERROR);
			return;
		}
		/* save obj mesh */
		mesh->saveAs(filename.ToStdString());

		/* If just points, need to correlate texcoords with vertices */
		if (mesh->n_faces() == 0) {
			FILE *out = fopen(filename.c_str(),"a");
			if (out) {
				for (unsigned int i = 1; i <= mesh->n_vertices(); i++) {
					fprintf(out,"f %d/%d\n",i,i);
				}
				fclose(out);
			}
		}
	}
}

void MainFrame::OnRefreshParWin(wxCommandEvent& WXUNUSED(event))
{
	parwin->partree->Reload();
}

void MainFrame::OnParProfile(wxCommandEvent& WXUNUSED(event))
{
	ParticleSystems *particlesystems = &(wxGetApp().particlesystems);
	wxString report = particlesystems->profile().c_str();

	wxFrame* frame=new wxFrame(NULL, -1, "Profile", wxDefaultPosition, wxSize(600,300));	
	wxTextCtrl * txtControl = new wxTextCtrl(frame, -1, "", wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE);
	txtControl->AppendText(report);
	frame->Show();
}

void MainFrame::OnPauseAll(wxCommandEvent& WXUNUSED(event))
{
	ParticleSystems *particlesystems = &(wxGetApp().particlesystems);
	particlesystems->pause();
}

void MainFrame::OnPlayAll(wxCommandEvent& WXUNUSED(event))
{
	ParticleSystems *particlesystems = &(wxGetApp().particlesystems);
	particlesystems->play();
}

void MainFrame::OnOpenGLExts(wxCommandEvent& WXUNUSED(event))
{
	wxString exts((char *)glGetString(GL_EXTENSIONS));
	exts.Replace(" ","\n");

#if defined(WIN32)
	typedef const char * (WINAPI * PFNWGLGETEXTENSIONSSTRINGARBPROC) (HDC hdc);
	PFNWGLGETEXTENSIONSSTRINGARBPROC wglGetExtensionsStringARB = 0;
	wglGetExtensionsStringARB = (PFNWGLGETEXTENSIONSSTRINGARBPROC)wglGetProcAddress("wglGetExtensionsStringARB");
	if(wglGetExtensionsStringARB)
	{
		wxString wgl_exts((char *)wglGetExtensionsStringARB(wglGetCurrentDC()));
		wgl_exts.Replace(" ","\n");
		exts += "\n" + wgl_exts;
	}
#endif

	//wxMessageBox(exts,"Available OpenGL Extensions",wxICON_INFORMATION | wxOK,this);

	// Use the log window to get a cheap textctrl to show all the 
	new wxLogWindow(this,"Available OpenGL Extensions",true,false);
	wxLogMessage(exts);
}

void MainFrame::OnZoomIn(wxCommandEvent& WXUNUSED(event))
{
	displaywin->zoom += 1.0;
	displaywin->SetupView(FALSE,0,0);
}

void MainFrame::OnZoomOut(wxCommandEvent& WXUNUSED(event))
{
	displaywin->zoom -= 1.0;
	displaywin->SetupView(FALSE,0,0);
}

void MainFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
	Close(TRUE);
	if(wxGetApp().logWindow) (wxGetApp().logWindow)->Destroy();

#ifdef WB_USE_CGAL
	if(wxGetApp().logWindow) (wxGetApp().clipArtFrame)->Destroy();
#endif

}

void MainFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
	wxMessageBox("sample code by John","About wxModeler",wxOK | wxICON_INFORMATION, this);
}

void MainFrame::OnIdle(wxIdleEvent &event)
{
	/* Redirect cerr into the string errors */
	if(wxGetApp().logWindow){
		cerrbuforig = std::cerr.rdbuf();
		std::cerr.rdbuf((wxGetApp().logWindow)->logTxtControl);
		coutbuforig = std::cout.rdbuf();
		std::cout.rdbuf((wxGetApp().logWindow)->logTxtControl);
	}
	wxGetApp().particlesystems.fullUpdate();
	displaywin->Refresh();
	event.RequestMore(true);	//ensure further OnIdles are called

	if(wxGetApp().logWindow){
		std::cerr.rdbuf(cerrbuforig);
		std::cout.rdbuf(coutbuforig);
	}
}

void MainFrame::OnSaveGL(wxCommandEvent& WXUNUSED(event))
{
	wxFileDialog *filedialog = new wxFileDialog(this, "Save to a File", "", "",
												"JPG (*.jpg)|*.jpg|"
												"Other files (*.*)|*.*", wxFD_SAVE );
	if (filedialog->ShowModal() == wxID_CANCEL)
		return;
	
	wxString filename(filedialog->GetPath());
	
	displaywin->CaptureGL(filename);
	
}

void MainFrame::OnSaveSVG(wxCommandEvent& WXUNUSED(event))
{
	
	wxFileDialog *filedialog = new wxFileDialog(this, "Save to a File", "", "",
												"SVG (*.svg)|*.svg|"
												"Other files (*.*)|*.*", wxFD_SAVE);
	if (filedialog->ShowModal() == wxID_CANCEL)
		return;
	
	wxString filename(filedialog->GetPath());
	
	int w, h;
    displaywin->GetClientSize(&w, &h);
	
	wxSVGFileDC svgDC (filename, w, h);
    
	displaywin->OnDraw(svgDC);
	
    
}

void MainFrame::OnStartSaveSVGAnimation(wxCommandEvent& WXUNUSED(event))
{
	wxGetApp().captureSVGflag = true;
}

void MainFrame::OnStopSaveSVGAnimation(wxCommandEvent& WXUNUSED(event))
{
	wxGetApp().captureSVGflag = false;
}


void MainFrame::OnStartSaveGLAnimation(wxCommandEvent& WXUNUSED(event))
{
	wxGetApp().captureCanvasFlag = true;
}

void MainFrame::OnStopSaveGLAnimation(wxCommandEvent& WXUNUSED(event))
{
	wxGetApp().captureCanvasFlag = false;
}


void MainFrame::OnShowLog(wxCommandEvent& WXUNUSED(event))
{
	//create the window log and the text control that captures the std::cout --ms
	if(wxGetApp().logWindow)
	{
		(wxGetApp().logWindow)->Show();
		return;
	}
	wxGetApp().logWindow = new LogFrame("Log Window", wxDefaultPosition, wxSize(600,300));
	(wxGetApp().logWindow)->Show();
}

void MainFrame::OnShowClipArt(wxCommandEvent& WXUNUSED(event))
{
#ifdef WB_USE_CGAL
	//create the clipart window  --ms
	if(wxGetApp().clipArtFrame)
	{
		(wxGetApp().clipArtFrame)->Show();
		return;
	}
	wxGetApp().clipArtFrame = new ClipArtFrame("Clip Art",wxDefaultPosition,wxSize(600,600));
	(wxGetApp().clipArtFrame)->Show();
#endif
}

