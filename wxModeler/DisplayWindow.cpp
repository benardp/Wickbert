#include "DisplayWindow.h"
#include "wxModeler.h"
#include "LogFrame.h"
#undef min
#undef max
#include "Particles/Attributes/ParticlePosition.h"
#include "Particles/ParticleShader.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include "wx/image.h"
#include <time.h>

#ifdef WB_USE_CGAL
#include "SVGshader.h"
#endif

BEGIN_EVENT_TABLE(DisplayWindow, wxGLCanvas)
	EVT_SIZE(DisplayWindow::OnSize)
	EVT_PAINT(DisplayWindow::OnPaint)
	EVT_MOUSE_EVENTS(DisplayWindow::OnMouseEvent)
	EVT_CHAR(DisplayWindow::OnKey)
	EVT_ERASE_BACKGROUND(DisplayWindow::OnEraseBackground)
//	EVT_KEY_DOWN(DisplayWindow::OnKeyDown)
END_EVENT_TABLE()

static int attribute_list[] = { WX_GL_RGBA, WX_GL_DOUBLEBUFFER, 0 };

DisplayWindow::DisplayWindow(wxWindow *parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style, const wxString& name)
    : wxGLCanvas(parent, id, attribute_list, pos, size, style|wxFULL_REPAINT_ON_RESIZE|wxSUNKEN_BORDER, name)
{
	selected_p = 0;
	selected_particle = -1;
	xrot = yrot = 0.0f;
	zoom = -5.0f;
	frame_count = 1;
	do_init = true; //initialize GL
	//I moved the GL initialization in initGL(), since 
	//wxGTK on Linux does not like initialization in the constructor of the canvas
	#ifdef WB_USE_CGAL
	svgShader = new SVGshader();
	#endif

    wxGLContextAttrs cxtAttrs;
    cxtAttrs.OGLVersion(2, 1).Robust().ResetIsolation().EndList();
    glContext = new wxGLContext(this, NULL, &cxtAttrs);
}

DisplayWindow::~DisplayWindow()
{
	#ifdef WB_USE_CGAL
	delete svgShader;
	#endif
}

void DisplayWindow::initGL()
{
	//SetCurrent();
	glClearColor(1.0,1.0,1.0,1.0);	
	glEnable(GL_DEPTH_TEST);
	static float ambient[] = {0.1f, 0.1f, 0.1f, 1.0f};
	static float diffuse[] = {1.0f, 1.0f, 1.0f, 1.0f};
	static float position0[] = {10.0f, 20.0f, 20.0f, 0.0f};
	static float position1[] = {0.0f, 0.0f, -20.0f, 0.0f};
	static float front_mat_shininess[] = {60.0f};
	static float front_mat_specular[] = {0.2f, 0.2f, 0.2f, 1.0f};
	static float front_mat_diffuse[] = {0.5f, 0.28f, 0.38f, 1.0f};	
	static float lmodel_ambient[] = {1.0f, 1.0f, 1.0f, 1.0f};
	static float lmodel_twoside[] = {GL_TRUE};	
	
	  /* speedups */
	glEnable(GL_DITHER);
	glShadeModel(GL_SMOOTH);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);

	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, position0);
	glEnable(GL_LIGHT0);
	
	glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT1, GL_POSITION, position1);
	glEnable(GL_LIGHT1);
	
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
	glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside);
	glEnable(GL_LIGHTING);
	
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, front_mat_shininess);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, front_mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, front_mat_diffuse);

	glEnable(GL_LINE_SMOOTH);	// antialias lines
}

void
DisplayWindow::SetupView(int picking, int x, int y)
{
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	int w,h;
    GetClientSize(&w, &h);

    if (!glContext) return;
    SetCurrent(*glContext);
	float aspect = (float)w/h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (picking) gluPickMatrix(x,viewport[3] - y,5.0,5.0,viewport);
	glFrustum(-0.5, 0.5, -0.5/aspect, 0.5/aspect, 0.5, 10.0 );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef( 0.0, 0.0, zoom);
	glRotatef(yrot,0.0,1.0,0.0);
	glRotatef(xrot,1.0,0.0,0.0);
}

/* This routine is just OpenGL, no WX (except wxGetApp) */
void DisplayWindow::render()
{
	ParticleSystems *ps = &(wxGetApp().particlesystems);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glInitNames();
	for (unsigned int j = 0; j < ps->size(); j++) {
		glPushName(j);
		for(unsigned int i = 0; i < (*ps)[j]->particles.size(); i++) {
			glPushName(i);
			(*ps)[j]->particles[i]->draw();
			glPopName();
		}
		glPopName();
	}

	glFlush();
	SwapBuffers();
	//added for GL canvas capturing -ms
	if(wxGetApp().captureCanvasFlag)
		this->captureCanvas();	
	ps->animation();
	if(wxGetApp().captureSVGflag)
		this->captureSVG();
}

void  DisplayWindow::CaptureGL(const wxString& filename)
{
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	int w, h;
	w = viewport[2];
	h = viewport[3];
	
	GLubyte* image = new GLubyte[w * h * 3];
	
	glPixelStorei(GL_PACK_ALIGNMENT, 1); // force 1 byte alignment
	
	glReadBuffer(GL_FRONT);
	
	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, image);
	
	///create an image with the content of the framebuffer
	wxImage::AddHandler( new wxJPEGHandler );
	wxImage wx_image(w,h,image);
	
	///save the image to a file
	wx_image.SaveFile(filename, wxBITMAP_TYPE_JPEG );
		
}
void  DisplayWindow::captureCanvas()
{

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	int w, h;
	w = viewport[2];
	h = viewport[3];
	
	GLubyte* image = new GLubyte[w * h * 3];
	
	glPixelStorei(GL_PACK_ALIGNMENT, 1); // force 1 byte alignment
	
	glReadBuffer(GL_FRONT);
	
	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, image);
	
	///create an image with the content of the framebuffer
	wxImage::AddHandler( new wxJPEGHandler );
	wxImage wx_image(w,h,image);
	
	///generate a new file name	
	wxString filename;
	wxString zeros;
	int digits = 0;
	while (int(frame_count / pow(10.0,digits)) != 0)
		digits++;
	for(int i = 0; i < 4 - digits; i++) zeros << 0;
	filename <<  zeros << frame_count << ".jpg";
	///save the image to a file
	wx_image.SaveFile(filename, wxBITMAP_TYPE_JPEG );
	frame_count++;
}

void DisplayWindow::captureSVG()
{
	///generate a new file name	
	wxString filename;
	wxString zeros;
	int digits = 0;
	while (int(frame_count / pow(10.0,digits)) != 0)
		digits++;
	for(int i = 0; i < 4 - digits; i++) zeros << 0;
	filename <<  zeros << frame_count << ".svg";
	///save the SVG to a file
	int w, h;
    this->GetClientSize(&w, &h);
	
	wxSVGFileDC svgDC (filename, w, h);
    
	this->OnDraw(svgDC);
   	frame_count++;
}

bool DisplayWindow::pick(int x, int y)
{
	GLuint selbuf[1024];
	glSelectBuffer(1024, selbuf);

	SetupView(TRUE,x,y);

	glRenderMode(GL_SELECT);
	render();

	SetupView(FALSE,0,0);
	GLuint hits = glRenderMode(GL_RENDER);
	if (hits > 1000000) hits = 0;
	/* selbuf now contains a number of records (= hits) of the form:
	 * ... 3 znear zfar ps-index particles-index particle-index ...
	 *
	 * need to find the nearest and select it
	 */

	GLuint closest = 0;
	int picked_ps = -1;
	int picked_p = -1;
	int picked_particle = -1;
	for (unsigned int i = 0; i < hits; i++) {
		if (selbuf[6*i + 1] < closest || i == 0) {
			closest = selbuf[6*i + 1];
			picked_ps = selbuf[6*i + 3];
			picked_p = selbuf[6*i + 4];
			picked_particle = selbuf[6*i + 5];
		}
	}

	selected_p = 0;
	selected_particle = -1;

	/* this keeps an empty pick() from resetting the current selection */
	if (!hits) return false;

	ParticleSystems *ps = &(wxGetApp().particlesystems);
	for (unsigned int j = 0; j < ps->size(); j++) {
		for (unsigned int i = 0; i < (*ps)[j]->particles.size(); i++) {
			if (j == picked_ps && i == picked_p) {
				(*ps)[j]->particles[i]->selectedParticle = picked_particle;
				selected_p = (*ps)[j]->particles[i];
				selected_particle = picked_particle;
			} else {
				(*ps)[j]->particles[i]->selectedParticle = -1;
			}
		}
	}
	return true;
}

void DisplayWindow::OnPaint( wxPaintEvent& event )
{
    // This is a dummy, to avoid an endless succession of paint messages.
    // OnPaint handlers must always create a wxPaintDC.
    wxPaintDC dc(this);

    if (!glContext) return;

    SetCurrent(*glContext);

// Initialize OpenGL
    if (do_init)
    {
        this->initGL();
        do_init = false;
    }

    render();
}

void DisplayWindow::OnSize(wxSizeEvent& event)
{
    // this is also necessary to update the context on some platforms
//    wxGLCanvas::OnSize(event);
	resize();
}

void DisplayWindow::resize()
{
    // set GL viewport (not called by wxGLCanvas::OnSize on all platforms...)
    int w, h;
    GetClientSize(&w, &h);
    if (!glContext) return;
    SetCurrent(*glContext);
    glViewport(0, 0, (GLint) w, (GLint) h);
	SetupView(FALSE,0,0);
}

void DisplayWindow::OnKey(wxKeyEvent& event)
{
	int key = event.GetKeyCode();
	switch (key) {
		case WXK_DELETE:
			if (selected_p && selected_particle != -1) {
				selected_p->removeParticle(selected_particle);
				selected_p = NULL; selected_particle = -1;
				wxGetApp().particlesystems.clearSelection();
			}
			break;
		default:
			event.Skip();
	}
}

void DisplayWindow::OnMouseEvent(wxMouseEvent& event)
{
	static int rot_dragging = 0;
	static int particle_dragging = 0;
	static float last_x, last_y;

	//printf("%f %f %d\n", event.GetX(), event.GetY(), (int)event.LeftIsDown());
	if (event.GetWheelRotation() != 0) {
		zoom += event.GetWheelRotation() / 300.0f;
		SetupView(FALSE,0,0);
		Refresh(TRUE);
	}
	if (event.MiddleIsDown())
	{
		zoom += (event.GetY() - last_y)/100;
		SetupView(FALSE,0,0);
		Refresh(TRUE);
	}
	else if(event.LeftIsDown()) {
		SetFocus();

		if (rot_dragging) {
			xrot += (event.GetY() - last_y)*1.0;
			yrot += (event.GetX() - last_x)*1.0;
			SetupView(FALSE,0,0);
			Refresh(TRUE);
		} else if (particle_dragging) {
			if (!selected_p || selected_particle == -1) {
				particle_dragging = 0;
				return;
			}

			float dx = event.GetX() - last_x;
			float dy = last_y - event.GetY();
			//pick(event.GetX(),event.GetY());

			GLdouble modelview[16];
			GLdouble projection[16];
			GLint viewport[4];
			glGetDoublev(GL_MODELVIEW_MATRIX,modelview);
			glGetDoublev(GL_PROJECTION_MATRIX,projection);
			glGetIntegerv(GL_VIEWPORT,viewport);

//			ParticlePosition *posattr = selected_p->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
			ParticlePosition *posattr;
			if (selected_p->findAttribute(posattr)) {
				gmVector3 pos = posattr->getPosition(selected_particle);

				GLdouble objx = pos[0], objy = pos[1], objz = pos[2];
				GLdouble winx,winy,winz;

				gluProject(objx,objy,objz,
							modelview,projection,viewport,
							&winx,&winy,&winz);

				winx += dx;
				winy += dy;

				gluUnProject(winx,winy,winz,
							modelview,projection,viewport,
							&objx,&objy,&objz);

				pos.assign(objx,objy,objz);
				posattr->setPosition(selected_particle, pos);
				//[tk 6.01.05] added this code to try to update the particle positions when I move a particle
				//[tk 6.01.05] seems to have done the trick - now particles moved with the mouse tell their ParticlePosition attribute that they have been updated
				// This is useful for RBF control particles.  I don't know why the regular Witkin-Heckbert control particles didn't require this kind of feedback.  Perhaps they are
				// updated like this at every frame regardless of movement.
				// It might be better to implement this the same way as the regular Witkin-Heckbert control particles - to update RBF centers the same way, I will have to find out how they are actually
				// implemented.
				posattr->changed = true;
			}

		} else {
			if (pick(event.GetX(),event.GetY()))
				particle_dragging = 1;
			else
				rot_dragging = 1;
		}
	} else if (event.LeftUp() && rot_dragging) {
		// clicked on a non-particle, clear selection
		rot_dragging = 0;
		selected_p = 0;
		selected_particle = -1;

		ParticleSystems *ps = &(wxGetApp().particlesystems);
		for (unsigned int j = 0; j < ps->size(); j++) {
			for (unsigned int i = 0; i < (*ps)[j]->particles.size(); i++) {
				(*ps)[j]->particles[i]->selectedParticle = -1;
			}
		}

	} else {
		rot_dragging = 0;
		particle_dragging = 0;
	}

	int e = 0;
	if (event.LeftDown()) {
		e = WB_CLICK;
	} else if (event.LeftDClick()) {
		e = WB_DOUBLECLICK;
	} else if (event.LeftUp()) {
		e = WB_DRAG;
	} else if (event.RightDown()) {
		e = WB_RIGHTCLICK;
	}
	if (event.ShiftDown()) e |= WB_SHIFT;
	if (event.ControlDown()) e |= WB_CTRL;
	if (event.AltDown()) e |= WB_ALT;

	ParticleSystems *ps = &(wxGetApp().particlesystems);

	if(wxGetApp().logWindow){
		(wxGetApp().mainframe)->cerrbuforig = std::cerr.rdbuf();
		std::cerr.rdbuf((wxGetApp().logWindow)->logTxtControl);
		(wxGetApp().mainframe)->coutbuforig = std::cout.rdbuf();
		std::cout.rdbuf((wxGetApp().logWindow)->logTxtControl);
	}

	for (unsigned int j = 0; j < ps->size(); j++){
		//pass the Euler angles and the zoom factor
		(*ps)[j]->setEulerAnglesAndZoom(xrot,yrot,zoom);
		for(unsigned int i = 0; i < (*ps)[j]->particles.size(); i++)
			(*ps)[j]->particles[i]->event(e);
	}

	if(wxGetApp().logWindow){
		std::cerr.rdbuf((wxGetApp().mainframe)->cerrbuforig);
		std::cout.rdbuf((wxGetApp().mainframe)->coutbuforig);
	}
	
	//update to new position
	last_x = event.GetX();
	last_y = event.GetY();

}

void DisplayWindow::OnEraseBackground(wxEraseEvent& event)
{
    // Do nothing, to avoid flashing.
}

void DisplayWindow::OnDraw(wxDC& dc)
{
#ifdef WB_USE_CGAL
	svgShader->setParticleSystems(&(wxGetApp().particlesystems));
	svgShader->drawSVG(dc); 
#endif
}

void DisplayWindow::OnDrawPlanarMap(wxDC& dc)
{
#ifdef WB_USE_CGAL
	svgShader->setParticleSystems(&(wxGetApp().particlesystems));
	clock_t start = clock();
	svgShader->builAndDrawPlanarMap(dc);
	clock_t finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;
	std::cout << "Planar Map duration (s): " << duration << std::endl;

#endif
}

void DisplayWindow::BuildSkeleton(wxDC& dc)
{
#ifdef WB_USE_CGAL
	svgShader->setParticleSystems(&(wxGetApp().particlesystems));
	svgShader->buildSkeleton(dc);
#endif
}

void  DisplayWindow::drawSVGnew(wxDC& dc)
{
#ifdef WB_USE_CGAL
	svgShader->setParticleSystems(&(wxGetApp().particlesystems));
	clock_t start = clock();
	svgShader->drawSVGnew(dc);	
	clock_t finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;
	std::cout << "Final drawing duration (s): " << duration << std::endl;
#endif
}

void  DisplayWindow::renderConstrainedDelaunay(wxDC& dc)
{
#ifdef WB_USE_CGAL
	svgShader->setParticleSystems(&(wxGetApp().particlesystems));
	clock_t start = clock();
	svgShader->buildAndDrawConstrainedDelaunay(dc);
	clock_t finish = clock();
	double duration = (double)(finish - start) / CLOCKS_PER_SEC;
	std::cout << "Constraind Delaunay duration (s): " << duration << std::endl;
#endif
}

/*
void  DisplayWindow::OnKeyDown(wxKeyEvent& event)
{
	Surfaces::iterator iter;
	for(iter = wxGetApp().surfaces.begin(); iter != wxGetApp().surfaces.end(); ++iter)
	{
		(*iter)->processEvent(event.GetKeyCode());
	}
	printf("%d /n",event.GetKeyCode());
	event.Skip();
}
*/
