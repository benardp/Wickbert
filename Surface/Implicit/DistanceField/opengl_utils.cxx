#include "opengl_utils.h"

const GLuint * nearest_hit(const GLuint buffer[], const GLint hits){
  assert(hits > 0);

  const GLuint * ptr = buffer;
  const GLuint * nearest = ptr;
  GLuint smallest_z1 = ptr[1];

  for(int i = 0; i < hits; i++){
    GLuint names = ptr[0];
    GLuint z1 = ptr[1];

    if(z1 < smallest_z1){
      nearest = ptr;
      smallest_z1 = z1;
    }

    ptr += 3 + names;
  }

  return nearest;
}



GLvoid set_material_color ( GLfloat r, GLfloat g, GLfloat b )
{
  GLfloat mat_specular[4] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mat_ambient_and_diffuse[4] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mat_shininess[1] = { 50.0 };
  
  mat_specular[0] = mat_ambient_and_diffuse[0] = r;
  mat_specular[1] = mat_ambient_and_diffuse[1] = g;
  mat_specular[2] = mat_ambient_and_diffuse[2] = b;

  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_ambient_and_diffuse);
}


scalar depth_z(int wx, int wy){
  GLfloat depth_z;
  GLint viewport[4];

  glGetIntegerv(GL_VIEWPORT, viewport);  
  glReadPixels(wx, viewport[3] - wy, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth_z);

  return depth_z;
}

Vector3d unproject_pixel(int wx, int wy, GLfloat depth_z){
	GLdouble modelMatrix[16];
  GLdouble projMatrix[16];
  GLint viewport[4];
  
  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
  glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
  glGetIntegerv(GL_VIEWPORT, viewport);
  
  // Notice that we have to correct the y pixel coordinate.  GL
  // assigns the origin to the lower left corner, while FLTK assigns
  // the origin to the upper left corner.
  
  GLdouble x,y,z;
  
  int status = gluUnProject(wx, viewport[3]-wy, depth_z,
			    modelMatrix, projMatrix, viewport,
			    &x, &y, &z);
  
  assert(status == GL_TRUE);
  
  return Vector3d((scalar) x, (scalar) y, (scalar) z);
}

Vector3d unproject_pixel(int wx, int wy)
{    
  return unproject_pixel(wx,wy,(GLfloat)depth_z(wx,wy));
}


void GLStatus::report_stacks()
{
    GLint depth;

    GLdouble modelMatrix[16];
    GLdouble projMatrix[16];
  
    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);


    glGetIntegerv(GL_PROJECTION_STACK_DEPTH, &depth);
    std::cerr << "   Projection stack depth = " << depth;
    glGetIntegerv(GL_MAX_PROJECTION_STACK_DEPTH, &depth);
    std::cerr << " (" << depth << " max)" << std::endl;

    for(int i = 0; i < 4; i++){
      for(int j = 0; j < 4; j++){
	std::cerr << projMatrix[i+4*j] << " ";
      }
      std::cerr << "\n";
    }
    std::cerr << "\n";

    glGetIntegerv(GL_MODELVIEW_STACK_DEPTH, &depth);
    std::cerr << "   ModelView stack depth = " << depth;
    glGetIntegerv(GL_MAX_MODELVIEW_STACK_DEPTH, &depth);
    std::cerr << " (" << depth << " max)" << std::endl;

    for(int i = 0; i < 4; i++){
      for(int j = 0; j < 4; j++){
	std::cerr << modelMatrix[i+4*j] << " ";
      }
      std::cerr << "\n";
    }
    std::cerr << "\n";


    glGetIntegerv(GL_TEXTURE_STACK_DEPTH, &depth);
    std::cerr << "   Texture stack depth = " << depth;
    glGetIntegerv(GL_MAX_TEXTURE_STACK_DEPTH, &depth);
    std::cerr << " (" << depth << " max)" << std::endl;
}

void GLStatus::report_errors(const char *msg)
{
    bool stack_error = false;

    for(GLenum err=glGetError(); err!=GL_NO_ERROR; err=glGetError())
    {
        std::cerr << "GL ERROR ";
        if( msg ) std::cerr << msg;
        std::cerr << ": " << (const char *)gluErrorString(err) << std::endl;

        if( err==GL_STACK_OVERFLOW || err==GL_STACK_UNDERFLOW )
            stack_error = true;
    }

    if( stack_error )  report_stacks();
}




void glTranslate(const Vector3d& v){
  glTranslatef((GLfloat)v.X(),(GLfloat)v.Y(),(GLfloat)v.Z());
}
void glVertex(const Vector3d& v){
  glVertex3f((GLfloat)v.X(),(GLfloat)v.Y(),(GLfloat)v.Z());
}
void glNormal(const Vector3d& n){
  glNormal3f((GLfloat)n.X(),(GLfloat)n.Y(),(GLfloat)n.Z());
}
