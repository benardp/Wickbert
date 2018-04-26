#ifndef VIEWWINDOW_H
#define VIEWWINDOW_H


#include "opengl.h"

#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/fl_file_chooser.H>
#include <FL/fl_ask.H>
#include <gfx/raster.h>

#include <arcballcontroller.h>

#include "selectable.h"
#include "triangle.h"

#include "adf.h"
#include "volumegraph.h"
#include "particle.h"

#include <set>
#include <vector>


typedef std::set< Selectable * > SelectableSet;

class DraggingController{
 private:
  SelectableSet selection;

  GLdouble modelMatrix[16];
  GLdouble projMatrix[16];
  GLint viewport[4];
  GLfloat depth_z;
 public:
  void mouse_down(int * where, const SelectableSet& s){
    assert(selection.empty());

    for(SelectableSet::const_iterator i = s.begin(); i != s.end(); i++){
      if((*i)->Selected()){
	selection.insert(*i);
      }
    }
    
    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
    glGetIntegerv(GL_VIEWPORT, viewport);
    
    glReadPixels(where[0], viewport[3] - where[1], 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth_z);
  }

  void mouse_drag(int *where, int *last){
    GLdouble x,y,z;
    
    int status;
    
    status = gluUnProject(last[0], viewport[3]-last[1], depth_z,
			  modelMatrix, projMatrix, viewport,
			  &x, &y, &z);
    Vector3d last_pos(x,y,z);
    assert(status == GL_TRUE);

    
    status = gluUnProject(where[0], viewport[3]-where[1], depth_z,
			  modelMatrix, projMatrix, viewport,
			  &x, &y, &z);
    Vector3d current_pos(x,y,z);
    assert(status == GL_TRUE);
    
    Vector3d delta = current_pos - last_pos;
    
    for(SelectableSet::iterator i = selection.begin(); i != selection.end(); i++){
      (*i)->Move(delta);
    }    
  }
  void mouse_up(){
    selection.clear();
  }
};


enum MouseMode{ SURFACE_SELECTION_MODE, BOX_SELECTION_MODE, CAMERA_MODE, DRAGGING_MODE, INACTIVE_MODE};

class DeformableMesh;
class LevelSet;

class ViewWindow : public Fl_Gl_Window {
public:
  int selected;

  SelectableSet selectables;
  
  ViewWindow(int x, int y, int w, int h);
  
  static ViewWindow * current;

  void Update();

  ADF * adf;

  bool snapshot_to_file(int format, std::string str_filename);

  DeformableMesh * deformable_mesh;


 private:
  ArcballController arcball;
  DraggingController dragging;
  
  //  DeformableMesh * deformable_mesh;

  void draw();
  void drawClear();
  void drawSetupTransform();
  void drawSetupLighting();
  void drawScene();
  void drawObjects();
  int handle(int e);

  MouseMode  BeginSurfaceSelection(int where[], int button);
  void      AdjustSurfaceSelection(int where[], int button);
  
};

#endif
