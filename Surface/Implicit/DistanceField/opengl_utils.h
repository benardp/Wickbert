#ifndef OPENGL_UTILS_H
#define OPENGL_UTILS_H

#include "vector3d.h"
#include "opengl.h"

const GLuint * nearest_hit(const GLuint buffer[], const GLint hits);
GLvoid set_material_color ( GLfloat r, GLfloat g, GLfloat b );
scalar depth_z(int wx, int wy);
Vector3d unproject_pixel(int wx, int wy);
Vector3d unproject_pixel(int wx, int wy, GLfloat depth_z);

class GLStatus{
  static void report_stacks();
  static void report_errors(const char *msg);
};


void glTranslate(const Vector3d&);
void glVertex(const Vector3d&);
void glNormal(const Vector3d&);


#endif
