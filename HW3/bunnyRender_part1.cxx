#include <stdio.h>
#include <iostream>
#include <math.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include "load_mesh.cpp"
#include "frame_timer.cpp"

void reshape(int w, int h);
void display();

int main(int argc, char** argv)
{
//glewInit();

// Initialize GLUT
glutInit(&argc,argv);
glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
glutInitWindowSize(512, 512);
glutInitWindowPosition(0,0);
glutCreateWindow("Bunny");
glutReshapeFunc(reshape);
// Initialize glew after window creation
glewInit();
init_timer();
glutDisplayFunc(display);

//Start rendering
glutMainLoop();
return 0;
}

void reshape(int w, int h)
  {
  glutPostRedisplay();
  }

void display()
  {
  start_timing();
  glEnable(GL_DEPTH_TEST);
  glLoadIdentity();
  gluLookAt(0,0,0,0,0,-1,0,1,0);
  glTranslatef(0.1, -1, -1.5);
  glScalef(10.0, 10.0, 10.0);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(-0.1, 0.1, -0.1, 0.1, 0.1, 1000);
  glViewport(0,0,512,512);

  glClearColor(0,0,0,0);
  glClearDepth(1000);
  glClear(GL_COLOR_BUFFER_BIT);
  glClear(GL_DEPTH_BUFFER_BIT);
  glDisable(GL_CULL_FACE);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);

  float lightColor[4] = {1, 1, 1, 0};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);
  float Ia[] = {0.2, 0.2, 0.2, 0};

  float sqrt3 = pow(3, 0.5);
  float p = 0;
  float l[] = {1/sqrt3, 1/sqrt3, 1/sqrt3, 0};
  float la[] = {0.0, 0.0, 0.0, 0};
  float ld[] = {1, 1, 1, 0};
  float ls[] = {0,0,0,0};
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Ia);
  glLightfv(GL_LIGHT0, GL_POSITION, l);
  glLightfv(GL_LIGHT0, GL_AMBIENT, la);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, ld);
  glLightfv(GL_LIGHT0, GL_SPECULAR, ls);
  glMaterialf(GL_FRONT, GL_SHININESS, p);

  float ka[] = {1, 1, 1, 0};
  float kd[] = {1, 1, 1, 0};
  float ks[] = {0,0,0,0};
  glMaterialfv(GL_FRONT, GL_AMBIENT, ka);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, kd);
  glMaterialfv(GL_FRONT, GL_SPECULAR, ks);

  int numVertices = 35947;
  int numTriangles = 69451;
  load_mesh("bunny.obj");
  glBegin(GL_TRIANGLES);
  for ( int i = 0; i < numTriangles; ++i )
    {
    int k0 = gTriangles[i].indices[0];
    int k1 = gTriangles[i].indices[1];
    int k2 = gTriangles[i].indices[2];
  
    glNormal3f(gNormals[k0].x,
               gNormals[k0].y,
               gNormals[k0].z);
    glVertex3f(gPositions[k0].x,
               gPositions[k0].y,
               gPositions[k0].z);

    glNormal3f(gNormals[k1].x,
               gNormals[k1].y,
               gNormals[k1].z);
    glVertex3f(gPositions[k1].x,
               gPositions[k1].y,
               gPositions[k1].z);

    glNormal3f(gNormals[k2].x,
               gNormals[k2].y,
               gNormals[k2].z);
    glVertex3f(gPositions[k2].x,
               gPositions[k2].y,
               gPositions[k2].z);
    }
  glEnd();

  float timeElapsed = stop_timing();
  gTotalFrames++;
  gTotalTimeElapsed += timeElapsed;
  float fps = gTotalFrames / gTotalTimeElapsed;
  char string[1024] = {0};
  sprintf(string, "OpenGL Bunny: %0.2f FPS", fps);
  std::cout << "Rendered in: " << fps << std::endl;
  glutSetWindowTitle(string);

  glutPostRedisplay();
  glutSwapBuffers();
  }
