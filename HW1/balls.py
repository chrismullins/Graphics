from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
from numpy import *
import numpy.linalg
from objectMath import *


name = 'BALLS'
ANTIALIASING_RAYS = 1
ANTIALIASING = 0
origin = array([0,0,0])
SQUARE = 262144
pixCenters = zeros((SQUARE, 3), dtype=float)
left = -0.1
right = 0.1
bottom = -0.1
top = 0.1
sphere1radius = 1
sphere1center = array([-4,0,-7])
sphere1ka = array([0.2,0.0,0.0])
sphere1kd = array([1.0,0.0,0.0])
sphere1ks = array([0.0,0.0,0.0])
sphere1p = 0
sphere2radius = 2
sphere2center = array([0,0,-7])
sphere2ka = array([0.0,0.2,0.0])
sphere2kd = array([0.0,0.5,0.0])
sphere2ks = array([0.5,0.5,0.5])
sphere2p = 32
sphere3radius = 1
sphere3center = array([4,0,-7])
sphere3ka = array([0.0,0.0,0.2])
sphere3kd = array([0.0,0.0,1.0])
sphere3ks = array([0.0,0.0,0.0])
sphere3p = 0
yplane = -2
planeka = array([0.2,0.2,0.2])
planekd = array([1.0,1.0,1.0])
planeks = array([0.0,0.0,0.0])
planepower = 0
lightPos = array([-4,4, -3])

def main():
	glutInit(sys.argv)
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH)
	glutInitWindowSize(512,512)
	glutCreateWindow(name)
	glutDisplayFunc(display)
	glMatrixMode(GL_MODELVIEW)
	glutMainLoop()
	return
def colorPixel(im, x, y, red, green, blue):
	DIM = 512
	index = y*DIM + x
	im[0][index*3] = red
	im[0][index*3+1] = green
	im[0][index*3+2] = blue

def rayTrace(i,j):
	intensity = array([0.0,0.0,0.0])
	u = pixCenters[j*512+i][0]
	v = pixCenters[j*512+i][1]
	w = pixCenters[j*512+i][2]
	if(ANTIALIASING == 1):
		u = u + ((random.random_sample() - 0.5)*0.5*0.2/512)
		v = v + ((random.random_sample() - 0.5)*0.5*0.2/512)
	pixCentersUpdate = array([u,v,w])
	intersectionSphere1 = sphereIntersect(origin, pixCentersUpdate, sphere1center, sphere1radius)
	intersectionSphere2 = sphereIntersect(origin, pixCentersUpdate, sphere2center, sphere2radius)
	intersectionSphere3 = sphereIntersect(origin, pixCentersUpdate, sphere3center, sphere3radius)
	planeIntersect = planeIntersection(origin, pixCentersUpdate, yplane)
	intersections = append(intersectionSphere1, intersectionSphere2)
	intersections = append(intersections, intersectionSphere3)
	intersections = append(intersections, planeIntersect)	
	if size(intersections) > 0:
		minValue = intersections.min()
		viewPlanePoint = array([pixCentersUpdate[0], pixCentersUpdate[1], pixCentersUpdate[2]])
		viewPlanePoint = viewPlanePoint / linalg.norm(viewPlanePoint)
		objectPoint = viewPlanePoint*minValue
		objectPointToLight = lightPos - objectPoint				
		if minValue in intersectionSphere1:
			intensity = computeSphereShading(u,v,minValue, sphere1center, lightPos, sphere1ka, sphere1kd, sphere1ks, sphere1p)
			sphere2LightIntersection = sphereIntersect(objectPoint, objectPointToLight, sphere2center, sphere2radius)
			sphere3LightIntersection = sphereIntersect(objectPoint, objectPointToLight, sphere3center, sphere3radius)
			spherePlaneIntersections = append(sphere2LightIntersection, sphere3LightIntersection)
			if size(spherePlaneIntersections) > 1:
				intensity = intensity / size(spherePlaneIntersections)
		elif minValue in intersectionSphere2:
			intensity = computeSphereShading(u,v,minValue, sphere2center, lightPos, sphere2ka, sphere2kd, sphere2ks, sphere2p)
			sphere1LightIntersection = sphereIntersect(objectPoint, objectPointToLight, sphere1center, sphere1radius)
			# this next line screws it up, I don't know why it's intersecting with the third sphere but it is
			#sphere3LightIntersection = sphereIntersect(objectPoint, objectPointToLight, sphere3center, sphere3radius)
			sphere3LightIntersection = array([])
			spherePlaneIntersections = append(sphere1LightIntersection, sphere3LightIntersection)
			if size(spherePlaneIntersections) > 1:
				intensity = intensity / size(spherePlaneIntersections)
		elif minValue in intersectionSphere3:
			intensity = computeSphereShading(u,v,minValue, sphere3center, lightPos, sphere3ka, sphere3kd, sphere3ks, sphere3p)
			sphere1LightIntersection = sphereIntersect(objectPoint, objectPointToLight, sphere1center, sphere1radius)
			sphere2LightIntersection = sphereIntersect(objectPoint, objectPointToLight, sphere2center, sphere2radius)
			spherePlaneIntersections = append(sphere1LightIntersection, sphere2LightIntersection)
			if size(spherePlaneIntersections) > 1:
				intensity = intensity / size(spherePlaneIntersections)
		elif minValue in planeIntersect:
			intensity = computePlaneShading(i,j,minValue,yplane,lightPos, planeka, planekd, planeks, planepower)
			sphere1LightIntersection = sphereIntersect(objectPoint, objectPointToLight, sphere1center, sphere1radius)
			sphere2LightIntersection = sphereIntersect(objectPoint, objectPointToLight, sphere2center, sphere2radius)
			sphere3LightIntersection = sphereIntersect(objectPoint, objectPointToLight, sphere3center, sphere3radius)
			spherePlaneIntersections = append(sphere1LightIntersection, sphere2LightIntersection)
			spherePlaneIntersections = append(spherePlaneIntersections, sphere3LightIntersection)
			if size(spherePlaneIntersections) > 1:
				#intensity = planeka*0xFF  # the way I did it below is a lot cooler, plus the shadows can intersect
				intensity = intensity / (size(spherePlaneIntersections))
			#colorPixel(image, i, j, intensity[0], intensity[1], intensity[2])
	return intensity	

def display():
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
	glutSwapBuffers()
	
	image = zeros((1, SQUARE*3), dtype=float)
	for i in range(0, SQUARE):
		image[0][i*3] = 0x00
		image[0][i*3+1] = 0x00
		image[0][i*3+2] = 0x00
	
	for i in range(0,512):
		for j in range(0,512):
			pixCenters[j*512+i][0] = left + ((right - left)*(i+0.5)/512.0)
			pixCenters[j*512+i][1] = bottom + ((top - bottom)*(j+0.5)/512.0)
			pixCenters[j*512+i][2] = -0.1
	for i in range(0, 512):
		for j in range(0, 512):
			if ANTIALIASING == 1:
				intensity = array([0.0,0.0,0.0])
				for k in range(0,ANTIALIASING_RAYS):
					intensity = intensity+rayTrace(i,j)
				intensity = intensity / ANTIALIASING_RAYS
				colorPixel(image, i, j, intensity[0], intensity[1], intensity[2])
			else:
				intensity = rayTrace(i,j)
				colorPixel(image, i, j, intensity[0], intensity[1], intensity[2])			

	image[0] = image[0] / 0xFF
	image[0] = image[0] ** (1/2.2)
	image[0] = image[0] * 0xFF
	glDrawPixels(512, 512, GL_RGB, GL_UNSIGNED_BYTE, image[0])
	glutSwapBuffers()
	return

if __name__=='__main__': main()






