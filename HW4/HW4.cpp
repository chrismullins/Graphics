#include <stdio.h>
#include <iostream>
#include <math.h>
#include <GL/glut.h>
#include <Eigen/Dense>
#include <cmath>


// declare data structures
float* pixels = NULL;
int ANTIALIASING = 1;
int RAYS_PER_PIXEL = 10;
int MAX_DISTANCE = 99999999;
int MAX_MIN_INTERSECTION = 9999999;
// for the intersections
const int HIT_NOTHING = 0;
const int HIT_RED_SPHERE = 1;
const int HIT_GREEN_SPHERE = 2;
const int HIT_BLUE_SPHERE = 3;
const int HIT_PLANE = 4;
static Eigen::Vector3f* origin = NULL;
int SQUARE = 262144;
static Eigen::Vector3f* pixCenters = NULL;
float left; float right; 
float bottom; float top;
float sphere1radius;
static Eigen::Vector3f* sphere1center;
static Eigen::Vector3f* sphere1ka;
static Eigen::Vector3f* sphere1kd;
static Eigen::Vector3f* sphere1ks;
float sphere1alpha;
float sphere1p;
float sphere2radius;
static Eigen::Vector3f* sphere2center;
static Eigen::Vector3f* sphere2ka;
static Eigen::Vector3f* sphere2kd;
static Eigen::Vector3f* sphere2ks;
float sphere2alpha;
float sphere2p;
float sphere3radius;
static Eigen::Vector3f* sphere3center;
static Eigen::Vector3f* sphere3ka;
static Eigen::Vector3f* sphere3kd;
static Eigen::Vector3f* sphere3ks;
float sphere3alpha;
float sphere3p;
float yplane;
static Eigen::Vector3f* planeka;
static Eigen::Vector3f* planekd;
static Eigen::Vector3f* planeks;
float planealpha;
float planepower = 0;
static Eigen::Vector3f* lightPos;

// function headers
void initialize_scene();
void colorPixel(int x, int y, float red, float green, float blue);
void display();
Eigen::Vector3f rayTrace(int i, int j);
Eigen::Vector3f intersectRayWithScene(Eigen::Vector3f origin, Eigen::Vector3f direction, int recursionDepth);
float intersectRayWithSphere(Eigen::Vector3f origin, 
                             Eigen::Vector3f direction,
                             Eigen::Vector3f sphereCenter,
                             float sphereRadius);
float intersectRayWithPlane(Eigen::Vector3f origin,
                            Eigen::Vector3f direction,
                            float yval);
Eigen::Vector3f computeSphereShading(Eigen::Vector3f origin,
                                     Eigen::Vector3f direction, 
                                     float intersectionDistance,
                                     Eigen::Vector3f center,
                                     float radius,
                                     Eigen::Vector3f ka,
                                     Eigen::Vector3f kd,
                                     Eigen::Vector3f ks,
                                     float p);
Eigen::Vector3f computePlaneShading(Eigen::Vector3f origin,
                                    Eigen::Vector3f direction,
                                    float intersectionDistance,
                                    float yval,
                                    Eigen::Vector3f ka,
                                    Eigen::Vector3f kd,
                                    Eigen::Vector3f ks,
                                    float p);

Eigen::Vector3f computePlaneShading(Eigen::Vector3f origin, Eigen::Vector3f direction,
                                    float intersectionDistance, float yval,
                                    Eigen::Vector3f ka, Eigen::Vector3f kd, Eigen::Vector3f ks, float p)
{
Eigen::Vector3f intensity;
Eigen::Vector3f normal;
normal << 0.0, 1.0, 0.0;
Eigen::Vector3f eyeToPlane = direction.normalized();
Eigen::Vector3f planePoint = intersectionDistance*eyeToPlane;
Eigen::Vector3f La = ka;
Eigen::Vector3f planePointToLight = lightPos[0] - planePoint;
planePointToLight.normalize();
planePointToLight.normalize();
Eigen::Vector3f Ld = kd*1*fmaxf(0, normal.dot(planePointToLight));
Eigen::Vector3f planeToEye = origin - planePoint;
planeToEye.normalize();
Eigen::Vector3f h = planeToEye + planePointToLight;
h.normalize();
Eigen::Vector3f Ls = ks*1*pow(fmaxf(0, normal.dot(h)),p);
intensity = La + Ld + Ls;
intensity(0) = fminf(1.0, intensity(0));
intensity(1) = fminf(1.0, intensity(1));
intensity(2) = fminf(1.0, intensity(2));
return intensity;
}

Eigen::Vector3f computeSphereShading(Eigen::Vector3f origin, Eigen::Vector3f direction, float intersectionDistance,
                                     Eigen::Vector3f center, float radius,
                                     Eigen::Vector3f ka, Eigen::Vector3f kd, Eigen::Vector3f ks, float p)
{
Eigen::Vector3f intensity;
Eigen::Vector3f La = ka;
Eigen::Vector3f vec = direction.normalized();
Eigen::Vector3f point = origin + intersectionDistance*vec;
Eigen::Vector3f normal = point - center;
normal.normalize();
Eigen::Vector3f l = lightPos[0] - point;
l.normalize();
Eigen::Vector3f Ld = kd*1*fmaxf(0, normal.dot(l));
Eigen::Vector3f v = origin - point;
v.normalize();
Eigen::Vector3f h = v+l;
h.normalize();
Eigen::Vector3f Ls = ks*1*pow(fmaxf(0,normal.dot(h)),p);
intensity = La + Ld + Ls;
intensity(0) = fminf(1.0, intensity(0));
intensity(1) = fminf(1.0, intensity(1));
intensity(2) = fminf(1.0, intensity(2));
return intensity;
}

float intersectRayWithPlane(Eigen::Vector3f origin, Eigen::Vector3f direction, float yval)
{
float tol = 0.001;
direction.normalize();
float far = 10000.0;
float near = 0.01;
float t = 1.0;
float hitDistance = 99999999;
Eigen::Vector3f proj = direction*far;
Eigen::Vector3f check = direction*t;
if (proj(1) > yval)
  {
  return hitDistance;
  }
else
  {
  // start the binary search
  while ( fabs(check(1) - yval) > tol)
    {
    if ( check(1) < yval )
      {
      far = t;
      t = near + ((far - near)/2);
      }
    else if ( check(1) > yval )
      {
      near = t;
      t = near +  ((far - near)/2);
      }
    check = direction*t;
    }
  hitDistance = t;
  return hitDistance;
  }
}

float intersectRayWithSphere(Eigen::Vector3f origin, Eigen::Vector3f direction, Eigen::Vector3f sphereCenter, float sphereRadius)
{
float minHitDistance = MAX_DISTANCE;
float tol = 0.0001;
direction.normalize();
Eigen::Vector3f dif = origin - sphereCenter;
float sqrtUnder = (pow(direction.dot(dif), 2)) - (direction.dot(direction) * dif.dot(dif) - pow(sphereRadius, 2));
if (sqrtUnder >= tol)
  {
  // 2 solutions exist
  float sol1 = (-(direction.dot(dif)) + (pow(sqrtUnder, 0.5))/(direction.dot(direction)));
  float sol2 = (-(direction.dot(dif)) - (pow(sqrtUnder, 0.5))/(direction.dot(direction)));
  minHitDistance = fminf(sol1, sol2);
  }
else if (sqrtUnder < tol && sqrtUnder > 0)
  {
  // 1 solution exists
  minHitDistance = (-(direction.dot(dif))/(direction.dot(direction)));
  }
if (minHitDistance < 0)
  {
  minHitDistance = MAX_DISTANCE;
  }
return minHitDistance;
}

Eigen::Vector3f intersectRayWithScene(Eigen::Vector3f origin, Eigen::Vector3f direction, int recursionDepth)
{
direction.normalize();
int intersectionEnum = HIT_NOTHING;
int numIntersections = 0;
float minIntersection = MAX_MIN_INTERSECTION;
float sphere1Intersection = intersectRayWithSphere(origin, direction, sphere1center[0], sphere1radius);
float alpha;
Eigen::Vector3f reflectedDirection;
Eigen::Vector3f hitNormal;
if (sphere1Intersection < minIntersection)
  {
  minIntersection = sphere1Intersection;
  intersectionEnum = 1;
  }
float sphere2Intersection = intersectRayWithSphere(origin, direction, sphere2center[0], sphere2radius);
if (sphere2Intersection < minIntersection)
  {
  minIntersection = sphere2Intersection;
  intersectionEnum = 2;
  }
float sphere3Intersection = intersectRayWithSphere(origin, direction, sphere3center[0], sphere3radius);
if (sphere3Intersection < minIntersection)
  {
  minIntersection = sphere3Intersection;
  intersectionEnum = 3;
  }
float planeIntersection = intersectRayWithPlane(origin, direction, yplane);
if (planeIntersection < minIntersection)
  {
  minIntersection = planeIntersection;
  intersectionEnum = 4;
  }
Eigen::Vector3f intensity;
intensity << 0.0, 0.0, 0.0;
Eigen::Vector3f hitOrigin = origin + minIntersection*direction;
if ( intersectionEnum != HIT_NOTHING)
  {
  Eigen::Vector3f hitToLight = lightPos[0] -  hitOrigin;
  hitToLight.normalize();
  if (intersectRayWithSphere(hitOrigin, hitToLight, sphere1center[0], sphere1radius) < MAX_DISTANCE && intersectionEnum != HIT_RED_SPHERE)
    {
    numIntersections++;
    }
  if (intersectRayWithSphere(hitOrigin, hitToLight, sphere2center[0], sphere2radius) < MAX_DISTANCE && intersectionEnum != HIT_GREEN_SPHERE)
    {
    numIntersections++;
    }
  if (intersectRayWithSphere(hitOrigin, hitToLight, sphere3center[0], sphere3radius) < MAX_DISTANCE && intersectionEnum != HIT_BLUE_SPHERE)
    {
    numIntersections++;
    }
  }
switch(intersectionEnum)
  {
  case HIT_NOTHING:
    break;
  case HIT_RED_SPHERE:
    intensity = computeSphereShading(origin, direction, minIntersection, sphere1center[0], sphere1radius, sphere1ka[0], sphere1kd[0], sphere1ks[0], sphere1p);
    alpha = sphere1alpha;
    hitNormal = hitOrigin - sphere1center[0];
    break;
  case HIT_GREEN_SPHERE:
    intensity = computeSphereShading(origin, direction, minIntersection, sphere2center[0], sphere2radius, sphere2ka[0], sphere2kd[0], sphere2ks[0], sphere2p);
    alpha = sphere2alpha;
    hitNormal = hitOrigin - sphere2center[0];
    break;
  case HIT_BLUE_SPHERE:
    intensity = computeSphereShading(origin, direction, minIntersection, sphere3center[0], sphere3radius, sphere3ka[0], sphere3kd[0], sphere3ks[0], sphere3p);
    alpha = sphere3alpha;
    hitNormal = hitOrigin - sphere3center[0];
    break;
  case HIT_PLANE:
    intensity = computePlaneShading(origin, direction, minIntersection, yplane, planeka[0], planekd[0], planeks[0], planepower);
    alpha = planealpha;
    hitNormal << 0.0, 1.0, 0.0;
    break;
  }
intensity = intensity / (numIntersections+1);
if (recursionDepth == 0)
  {
  return intensity;
  }
else
  {
  if (intersectionEnum == 0)
    {
    intensity << 0.0, 0.0, 0.0;
    return intensity;
    }
  hitNormal.normalize();
  Eigen::Vector3f v = origin - hitOrigin;
  v.normalize();
  reflectedDirection = 2*(hitNormal.dot(v))*hitNormal - v;
  reflectedDirection.normalize();
  return (1 - alpha)*intensity + alpha*intersectRayWithScene(hitOrigin, reflectedDirection, recursionDepth - 1);
  }
}

Eigen::Vector3f rayTrace(int i, int j)
{
Eigen::Vector3f intensity;
Eigen::Vector3f eye;
intensity << 0.0, 0.0, 0.0;
eye << 0.0, 0.0, 0.0;
float u = pixCenters[j*512+i](0);
float v = pixCenters[j*512+i](1);
float w = pixCenters[j*512+i](2);
Eigen::Vector3f pixCentersUpdate;
pixCentersUpdate << u, v, w;
Eigen::Vector3f direction = pixCentersUpdate.normalized();
if (ANTIALIASING == 1)
  {
  for(int i = 0; i < RAYS_PER_PIXEL; ++i)
    {
    Eigen::Vector3f rayDirection;
    float r = ((float) rand() / (RAND_MAX));
    rayDirection << u + (r - 0.5)*0.5*0.2/512,
                    v + (r - 0.5)*0.5*0.2/512,
                    w;
    rayDirection.normalize();
    intensity = intensity + intersectRayWithScene(origin[0], rayDirection, 2);
    }
  intensity = intensity / RAYS_PER_PIXEL;
  }
else
  {
  intensity = intersectRayWithScene(origin[0], direction, 2);
  }
return intensity;
}

void display()
{
glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
glutSwapBuffers();
for(int i = 0; i < 512; ++i)
  {
  for(int j = 0; j < 512; ++j)
    {
    pixCenters[j*512+i](0) = left + ((right - left)*(i+0.5)/512.0);
    pixCenters[j*512+i](1) = bottom + ((top - bottom)*(j+0.5)/512.0);
    pixCenters[j*512+i](2) = -0.1;
    }
  }
#pragma omp parallel for
for(int i = 0; i < 512; ++i)
  {
  for (int j = 0; j < 512; ++j)
    {
    Eigen::Vector3f intensity = rayTrace(i,j);
    for (int k = 0; k < 3; ++k)
      {
      intensity[k] = pow(intensity[k], 1/2.2);
      }
    colorPixel(i, j, intensity[0], intensity[1], intensity[2]);
    }
  }
glDrawPixels(512, 512, GL_RGB, GL_FLOAT, pixels);
glutSwapBuffers();
return;
}

void colorPixel(int x, int y, float red, float green, float blue)
{
int index = y*512 + x;
pixels[index*3] = red;
pixels[index*3+1] = green;
pixels[index*3+2] = blue;
}

void initialize_scene()
{
pixels = new float[512*512*3];
for (int i = 0; i < SQUARE*3; ++i)
  {
  pixels[i] = 0.0;
  }
origin = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
origin[0] << 0.0, 0.0, 0.0;

pixCenters = (Eigen::Vector3f *) malloc(SQUARE * sizeof(Eigen::Vector3f));
left = -0.1;
right = 0.1;
bottom = -0.1;
top = 0.1;

sphere1radius = 1.0;
sphere1center = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
sphere1ka = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
sphere1kd = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
sphere1ks = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
sphere1center[0] << -4, 0, -7;
sphere1ka[0] << 0.2, 0.0, 0.0;
sphere1kd[0] << 1.0, 0.0, 0.0;
sphere1ks[0] << 0.0, 0.0, 0.0;
sphere1alpha = 0.0;
sphere1p = 0.0;

sphere2radius = 2.0;
sphere2center = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
sphere2ka = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
sphere2kd = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
sphere2ks = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
sphere2center[0] << 0.0, 0.0, -7;
sphere2ka[0] << 0.0, 0.2, 0.0;
sphere2kd[0] << 0.0, 0.5, 0.0;
sphere2ks[0] << 0.5, 0.5, 0.5;
sphere2alpha = 0.0;
sphere2p = 32.0;

sphere3radius = 1.0;
sphere3center = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
sphere3ka = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
sphere3kd = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
sphere3ks = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
sphere3center[0] << 4.0, 0.0, -7.0;
sphere3ka[0] << 0.0, 0.0, 0.2;
sphere3kd[0] << 0.0, 0.0, 1.0;
sphere3ks[0] << 0.0, 0.0, 0.0;
sphere3alpha = 0.8;
sphere3p = 0.0;

planeka = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
planekd = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
planeks = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
planeka[0] << 0.2, 0.2, 0.2;
planekd[0] << 1.0, 1.0, 1.0;
planeks[0] << 0.0, 0.0, 0.0;
planealpha = 0.5;
planepower = 0.0;

lightPos = (Eigen::Vector3f *) malloc(sizeof(Eigen::Vector3f));
lightPos[0] << -4.0, 4.0, -3.0;

yplane = -2.0;
}

int main(int argc, char** argv)
{
glutInit(&argc, argv); 
glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
glutInitWindowSize(512,512);
glutInitWindowPosition(0,0); 
glutCreateWindow("balls PA4"); 
initialize_scene();
glutDisplayFunc(display);
glMatrixMode(GL_MODELVIEW);
glutMainLoop();
return 0;
}
