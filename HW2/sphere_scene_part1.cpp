//
//  sphere_scene.c
//  Rasterizer
//
//

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <GL/glut.h>
#include <Eigen/Dense>

int     gNumVertices    = 0;    // Number of 3D vertices.
int     gNumTriangles   = 0;    // Number of triangles.
int*    gIndexBuffer    = NULL; // Vertex indices for the triangles.
static Eigen::Vector3f* gVertices = NULL;
static Eigen::Vector4f* modelTransformed = NULL;
static Eigen::Vector4f* perspectiveTransformed = NULL;
static Eigen::Vector4f* viewportTransformed = NULL;
float* pixels = NULL;
float* depthBuffer = NULL;

Eigen::Matrix3f triangleViewportCoords(int triangleIndex);
float beta(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c, int x, int y);
float betax(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c, int x, int y);
float betay(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c, int x, int y);
float gamma(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c, int x, int y);
float gammax(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c, int x, int y);
float gammay(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c, int x, int y);
void colorPixel(int x, int y, float red, float green, float blue);

void display(float pixels[])
  {
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glutSwapBuffers();
  
  glDrawPixels(512, 512, GL_RGB, GL_FLOAT, pixels);
  glutSwapBuffers();
  return;
  }

void colorTriangle(int i)
  {
  Eigen::Matrix3f tVertices = triangleViewportCoords(i);
  Eigen::Vector3f v1; 
  Eigen::Vector3f v2;
  Eigen::Vector3f v3;
  v1 << tVertices(0,0), tVertices(0,1), tVertices(0,2);
  v2 << tVertices(1,0), tVertices(1,1), tVertices(1,2);
  v3 << tVertices(2,0), tVertices(2,1), tVertices(2,2);
  int xmin = floor(fminf(v1(0), fminf(v2(0), v3(0))));
  int xmax = ceil(fmaxf(v1(0), fmaxf(v2(0), v3(0))));
  int ymin = floor(fminf(v1(1), fminf(v2(1), v3(1))));
  int ymax = ceil(fmaxf(v1(1), fmaxf(v2(1), v3(1))));
  float b = beta(v1, v2, v3, xmin, ymin);
  float bx = beta(v1,v2,v3,1,0) - beta(v1,v2,v3,0,0);
  float by = beta(v1,v2,v3,0,1) - beta(v1,v2,v3,0,0);
  float g = gamma(v1, v2, v3,xmin,ymin);
  float gx = gamma(v1,v2,v3,1,0) - gamma(v1,v2,v3,0,0);
  float gy = gamma(v1,v2,v3,0,1) - gamma(v1,v2,v3,0,0);
  int n = xmax - xmin;
  for (int y = ymin; y <= ymax; ++y)
    {
    for (int x = xmin; x <=xmax; ++x)
      {
      b = beta(v1,v2,v3,x,y);
      g = gamma(v1,v2,v3,x,y);
      if (b >= 0 && g >= 0 && b+g <= 1)
        {
        colorPixel(x, y, 0.0, 1.0, 0.0);
        }
      }
    }
  return;
  }

float beta(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c, int x, int y)
  {
  float xa = a(0);
  float xb = b(0);
  float xc = c(0);
  float ya = a(1);
  float yb = b(1);
  float yc = c(1);
  float numerator = (ya - yc)*x + (xc - xa)*y + (xa*yc) - (xc*ya);
  float denominator = (ya - yc)*xb + (xc - xa)*yb + (xa*yc) - (xc*ya);
  float beta = numerator / denominator;
  return beta;
  }

float betax(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c, int x, int y)
  {
  return beta(a,b,c,x+1,y) - beta(a,b,c,x,y);
  }

float betay(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c, int x, int y)
  {
  return beta(a,b,c,x,y+1) - beta(a,b,c,x,y);
  }

float gamma(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c, int x, int y)
  {
  float xa = a(0);
  float xb = b(0);
  float xc = c(0);
  float ya = a(1);
  float yb = b(1);
  float yc = c(1);
  float numerator = (ya - yb)*x + (xb - xa)*y + xa*yb - xb*ya;
  float denominator = (ya - yb)*xc + (xb - xa)*yc + xa*yb - xb*ya;
  return numerator / denominator;
  }

float gammax(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c, int x, int y)
  {
  return gamma(a,b,c,x+1,y) - gamma(a,b,c,x,y);
  }

float gammay(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c, int x, int y)
  {
  return gamma(a,b,c,x,y+1) - gamma(a,b,c,x,y);
  }

Eigen::Matrix3f triangleVertices(int triangleIndex)
  {
  Eigen::Matrix3f vertices;
  int k0 = gIndexBuffer[3*triangleIndex + 0];
  int k1 = gIndexBuffer[3*triangleIndex + 1];
  int k2 = gIndexBuffer[3*triangleIndex + 2];
  Eigen::Vector3f v1;
  Eigen::Vector3f v2;
  Eigen::Vector3f v3;
  v1 << gVertices[k0](0), gVertices[k0](1), gVertices[k0](2);
  v2 << gVertices[k1](0), gVertices[k1](1), gVertices[k1](2);
  v3 << gVertices[k2](0), gVertices[k2](1), gVertices[k2](2);
  vertices << v1(0),v1(1),v1(2),v2(0),v2(1),v2(2),v3(0),v3(1),v3(2);
  return vertices;
  }

Eigen::Matrix3f triangleWorldCoords(int triangleIndex)
  {
  Eigen::Matrix3f vertices;
  int k0 = gIndexBuffer[3*triangleIndex + 0];
  int k1 = gIndexBuffer[3*triangleIndex + 1];
  int k2 = gIndexBuffer[3*triangleIndex + 2];
  Eigen::Vector3f v1;
  Eigen::Vector3f v2;
  Eigen::Vector3f v3;
  v1 << modelTransformed[k0](0), modelTransformed[k0](1), modelTransformed[k0](2);
  v2 << modelTransformed[k1](0), modelTransformed[k1](1), modelTransformed[k1](2);
  v3 << modelTransformed[k2](0), modelTransformed[k2](1), modelTransformed[k2](2);
  vertices << v1(0),v1(1),v1(2),v2(0),v2(1),v2(2),v3(0),v3(1),v3(2);
  return vertices;
  }

Eigen::Matrix3f triangleViewportCoords(int triangleIndex)
  {
  Eigen::Matrix3f vertices;
  int k0 = gIndexBuffer[3*triangleIndex + 0];
  int k1 = gIndexBuffer[3*triangleIndex + 1];
  int k2 = gIndexBuffer[3*triangleIndex + 2];
  Eigen::Vector3f v1;
  Eigen::Vector3f v2;
  Eigen::Vector3f v3;
  v1 << viewportTransformed[k0](0), viewportTransformed[k0](1), viewportTransformed[k0](2);
  v2 << viewportTransformed[k1](0), viewportTransformed[k1](1), viewportTransformed[k1](2);
  v3 << viewportTransformed[k2](0), viewportTransformed[k2](1), viewportTransformed[k2](2);
  vertices << v1(0),v1(1),v1(2),v2(0),v2(1),v2(2),v3(0),v3(1),v3(2);
  return vertices;
  }

void colorPixel(int x, int y, float red, float green, float blue)
  {
  int DIM = 512;
  int index = 512*y + x;
  pixels[index*3] = red;
  pixels[index*3+1] = green;
  pixels[index*3+2] = blue;
  }

void create_scene()
{
    int width   = 32;
    int height  = 16;
    
    float theta, phi;
    int t;
    
    gNumVertices    = (height - 2) * width + 2;
    gNumTriangles   = (height - 2) * (width - 1) * 2;
    
    gVertices = (Eigen::Vector3f *) malloc(gNumVertices * sizeof(Eigen::Vector3f));

    gIndexBuffer    = new int[3*gNumTriangles];
    
    t = 0;
    for (int j = 1; j < height-1; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            theta = (float) j / (height-1) * M_PI;
            phi   = (float) i / (width-1)  * M_PI * 2;
            
            float   x   = sinf(theta) * cosf(phi);
            float   y   = cosf(theta);
            float   z   = -sinf(theta) * sinf(phi);
            
            gVertices [t](0) = x;
            gVertices [t](1) = y;
            gVertices [t](2) = z;
            
            t++;
        }
    }
    
    gVertices [t](0) = 0;
    gVertices [t](1) = 1;
    gVertices [t](2) = 0;    

    t++;
    
    gVertices [t][0] = 0;
    gVertices [t](1) = -1;
    gVertices [t](2) = 0;    

    t++;
    
    t = 0;
    for (int j = 0; j < height-3; ++j)
    {
        for (int i = 0; i < width-1; ++i)
        {
            gIndexBuffer[t++] = j*width + i;
            gIndexBuffer[t++] = (j+1)*width + (i+1);
            gIndexBuffer[t++] = j*width + (i+1);
            gIndexBuffer[t++] = j*width + i;
            gIndexBuffer[t++] = (j+1)*width + i;
            gIndexBuffer[t++] = (j+1)*width + (i+1);
        }
    }
    for (int i = 0; i < width-1; ++i)
    {
        gIndexBuffer[t++] = (height-2)*width;
        gIndexBuffer[t++] = i;
        gIndexBuffer[t++] = i + 1;
        gIndexBuffer[t++] = (height-2)*width + 1;
        gIndexBuffer[t++] = (height-3)*width + (i+1);
        gIndexBuffer[t++] = (height-3)*width + i;
    }
    
    Eigen::Matrix4f modelTranslate;
    modelTranslate << 1,0,0,0,
                      0,1,0,0,
                      0,0,1,-7,
                      0,0,0,1;
    Eigen::Matrix4f modelScale;
    modelScale << 2,0,0,0,
                  0,2,0,0,
                  0,0,2,0,
                  0,0,0,1;
   Eigen::Matrix4f modelTransform = modelTranslate*modelScale;
   Eigen::MatrixXd eVertices(4,gNumVertices);
   for (int i = 0; i < gNumVertices; ++i)
      {
      eVertices(0,i) = gVertices[i](0);
      eVertices(1,i) = gVertices[i](1);
      eVertices(2,i) = gVertices[i](2);
      eVertices(3,i) = 1;
      }

  modelTransformed = (Eigen::Vector4f *) malloc(gNumVertices * sizeof(Eigen::Vector4f));
  for (int i = 0; i < gNumVertices; ++i)
    {
    Eigen::Vector4f v;
    v << eVertices(0,i), eVertices(1,i), eVertices(2,i), eVertices(3,i);
    Eigen::Vector4f res;
    res = modelTransform*v;
    modelTransformed[i](0) = res(0);
    modelTransformed[i](1) = res(1);
    modelTransformed[i](2) = res(2);
    modelTransformed[i](3) = res(3);
    }

  float r = 0.1;
  float l = -0.1;
  float b = -0.1;
  float top = 0.1;
  float n = -0.1;
  float f = -1000.0;
  Eigen::Matrix4f perspectiveTransform;
  perspectiveTransform << ((2*n)/(r-l)), 0, (l+r)/(l-r), 0,
                          0, ((2*n)/(top-b)), (b+top)/(b-top), 0,
                          0, 0, (f+n)/(n-f), ((2*f*n)/(f-n)),
                          0,0,1,0;
  perspectiveTransformed = (Eigen::Vector4f *) malloc(gNumVertices * sizeof(Eigen::Vector4f));
  for (int i = 0; i < gNumVertices; ++i)
    {
    Eigen::Vector4f v;
    v << modelTransformed[i](0), modelTransformed[i](1), modelTransformed[i](2), modelTransformed[i](3);
    Eigen::Vector4f res;
    res = perspectiveTransform*v;
    perspectiveTransformed[i](0) = res(0);
    perspectiveTransformed[i](1) = res(1);
    perspectiveTransformed[i](2) = res(2);
    perspectiveTransformed[i](3) = res(3);
    }

  int nx = 512;
  int ny = 512;
  Eigen::Matrix4f viewportTransform;
  viewportTransform << (nx/2), 0, 0, (nx-1)/2,
                       0, (ny/2), 0, (ny-1)/2,
                       0, 0, 1, 0,
                       0, 0, 0, 1;

   viewportTransformed = (Eigen::Vector4f *) malloc(gNumVertices * sizeof(Eigen::Vector4f));
  for ( int i = 0; i < gNumVertices; ++i)
    {
    Eigen::Vector4f v;
    v << perspectiveTransformed[i](0), perspectiveTransformed[i](1), perspectiveTransformed[i](2), perspectiveTransformed[i](3);
    Eigen::Vector4f res;
    res = viewportTransform*v;
    float w = res(3);
    viewportTransformed[i](0) = res(0)/w;
    viewportTransformed[i](1) = res(1)/w;
    viewportTransformed[i](2) = res(2)/w;
    viewportTransformed[i](3) = res(3)/w;
    }
 
  Eigen::Matrix4f finalTransform = viewportTransform * perspectiveTransform * modelTransform;
  for( int i = 0; i < gNumVertices; ++i)
    {
    Eigen::Vector4f v;
    v << gVertices[i](0), gVertices[i](1), gVertices[i](2), 1.0;
    Eigen::Vector4f res = finalTransform*v;
    float w = res(3);
    viewportTransformed[i](0) = res(0)/w;
    viewportTransformed[i](1) = res(1)/w;
    viewportTransformed[i](2) = res(2)/w;
    viewportTransformed[i](3) = res(3)/w;
    }

  int size = 512*512;
  depthBuffer = new float[size];
  pixels = new float[size*3];
  for (int i = 0; i < size*3; ++i)
    {
    pixels[i] = 0.0;
    }
  for (int i = 0; i < gNumTriangles; ++i)
    {
     colorTriangle(i);
    }
  std::cout << "CHECK3" << std::endl;
  display(pixels);
  
}

int main(int argc, char** argv)
  {
  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(512,512);
  glutInitWindowPosition(0,0); 
  glutCreateWindow("No shading"); 
  create_scene();
  glMatrixMode(GL_MODELVIEW);
  glutMainLoop();

  return 0;
  }
