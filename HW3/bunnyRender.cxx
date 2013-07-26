#include <stdio.h>
#include <iostream>
#include <math.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include "load_mesh.cpp"
#include "frame_timer.cpp"

int numVertices = 35947;
int numTriangles = 69451;

void reshape(int w, int h);
void display();
void bufferData();

int main(int argc, char** argv)
{
//glewInit(); // Don't do this here, wait until after the window is initialized

// Initialize GLUT
glutInit(&argc,argv);
glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
glutInitWindowSize(512, 512);
glutInitWindowPosition(0,0);
glutCreateWindow("Bunny");
glutReshapeFunc(reshape);
// Initialize glew after window creation
glewExperimental=GL_TRUE;
glewInit();
init_timer();
bufferData();
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
  float l[] = {1, 1, 1, 0};
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

  
  glDrawElements(GL_TRIANGLES, 3*numTriangles+1, GL_UNSIGNED_INT, 0);
  // Don't uncomment these last two lines
  // Need to find a better place to do this
  //glDisableClientState(GL_VERTEX_ARRAY);
  //glDisableClientState(GL_NORMAL_ARRAY);

  float timeElapsed = stop_timing();
  gTotalFrames++;
  gTotalTimeElapsed += timeElapsed;
  float fps = gTotalFrames / gTotalTimeElapsed;
  char string[1024] = {0};
  sprintf(string, "OpenGL Bunny: %0.2f FPS", fps);
  //std::cout << "Rendered in: " << fps << std::endl;
  glutSetWindowTitle(string);

  glutPostRedisplay();
  glutSwapBuffers();
  }

void bufferData()
  {
  
  load_mesh("bunny.obj");

  GLfloat vertices[3*numVertices];
  GLfloat normals[3*numVertices];
  GLuint indices[3*numTriangles];
  for (int i = 0; i < numVertices; ++i)
    {
    vertices[3*i + 0] = gPositions[i].x;
    vertices[3*i + 1] = gPositions[i].y;
    vertices[3*i + 2] = gPositions[i].z;
    normals[3*i + 0] = gNormals[i].x;
    normals[3*i + 1] = gNormals[i].y;
    normals[3*i + 2] = gNormals[i].z;
    }
  std::cout << "Init positions and normals" << std::endl;
  for (int i = 0; i < numTriangles; ++i)
    {
    int k0 = gTriangles[i].indices[0];
    int k1 = gTriangles[i].indices[1];
    int k2 = gTriangles[i].indices[2];
    indices[3*i + 0] = k0;
    indices[3*i + 1] = k1;
    indices[3*i + 2] = k2;
    }
  std::cout << "Filled indices buffer" << std::endl;

  // This will create an interleaved array
  /*
  GLfloat data[3*2*numVertices];
  for (int i = 0; i < numVertices; ++i)
    {
    data[6*i + 0] = vertices[3*i + 0];
    data[6*i + 1] = vertices[3*i + 1];
    data[6*i + 2] = vertices[3*i + 2];
    data[6*i + 3] = normals[3*i + 0];
    data[6*i + 4] = normals[3*i + 1];
    data[6*i + 5] = normals[3*i + 2];
    }
  std::cout << "Filled data buffer" << std::endl;
  */   

  GLuint gVAO;

  GLuint gVBO[2];
  glGenBuffers(2,gVBO);
  GLuint positionBufferHandle = gVBO[0];
  GLuint normalBufferHandle = gVBO[1];

  glBindBuffer(GL_ARRAY_BUFFER, positionBufferHandle);
  /*
   *           GPU buffer you want to bind to
   *           |                Amount of data you want to put there
   *           |                |                              Name of local data aray
   *           |                |                              |         We don't want to change it
   *           |                |                              |         |                        */
  glBufferData(GL_ARRAY_BUFFER, 3*numVertices*sizeof(GLfloat), vertices, GL_STATIC_DRAW);

  glBindBuffer(GL_ARRAY_BUFFER, normalBufferHandle);
  glBufferData(GL_ARRAY_BUFFER, 3*numVertices*sizeof(GLfloat), normals, GL_STATIC_DRAW);

  glGenVertexArrays(1, &gVAO);
  glBindVertexArray(gVAO);
  std::cout << "Enabled VAO" << std::endl;

  /*
   * Replace glEnableVertexAttribArray with glEnableClientState
   * to work on Nvidia hardware
   */
  //glEnableVertexAttribArray(0);
  glEnableClientState(GL_VERTEX_ARRAY);
  //glEnableVertexAttribArray(1);
  glEnableClientState(GL_NORMAL_ARRAY);
 
  /*
   * Replace calls to glVertexAttribPointer with
   * gl<Vertex, Normal>Pointer
   */
  glBindBuffer(GL_ARRAY_BUFFER,positionBufferHandle);
  //glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,0);
  glVertexPointer(3, GL_FLOAT, 0, 0);

  glBindBuffer(GL_ARRAY_BUFFER,normalBufferHandle);
  //glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,0);
  glNormalPointer(GL_FLOAT, 0, 0);

  GLuint gEBO;
  glGenBuffers(1, &gEBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gEBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*numTriangles*sizeof(GLuint), indices, GL_STATIC_DRAW);
  std::cout << "Enabled EBO and buffered data with error code of " << glGetError() << std::endl;

  glBindVertexArray(gVAO);
  
  }
