// ==========================================================================
// Barebones OpenGL Core Profile Boilerplate
//    using the GLFW windowing system (http://www.glfw.org)
//
// Loosely based on
//  - Chris Wellons' example (https://github.com/skeeto/opengl-demo) and
//  - Camilla Berglund's example (http://www.glfw.org/docs/latest/quick.html)
//
// Author:  Sonny Chan, University of Calgary
// Date:    December 2015
// ==========================================================================

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <iterator>
#include <glm/glm.hpp>
#include "ImageWriter.hpp"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// Specify that we want the OpenGL core profile before including GLFW headers
#ifndef LAB_LINUX
#include <glad/glad.h>
#else
#define GLFW_INCLUDE_GLCOREARB
#define GL_GLEXT_PROTOTYPES
#endif
#include <GLFW/glfw3.h>

using namespace std;

// --------------------------------------------------------------------------
// constants and global vars

const float pi = 3.14159265359;
const int maxShapes = 32;
const int wWidth = 1024;
const int wHeight = 1024;
const bool antialiasing = false;

float fovDegrees = 38.0;
float fov = fovDegrees * pi / 180.0;
float diffRatio = 0.5;

int currScene = 2;
bool initialized = false;

glm::vec3 origin = glm::vec3(-0.01, 0.01, 0.0);
glm::vec3 viewRays[wWidth*wHeight];

GLfloat vertices[wWidth*wHeight][2];
GLfloat colours[wWidth*wHeight][3];

// --------------------------------------------------------------------------
// OpenGL utility and support function prototypes

void QueryGLVersion();
bool CheckGLErrors();

string LoadSource(const string &filename);
GLuint CompileShader(GLenum shaderType, const string &source);
GLuint LinkProgram(GLuint vertexShader, GLuint fragmentShader);

// --------------------------------------------------------------------------
// Functions to set up OpenGL shader programs for rendering

struct MyShader
{
	// OpenGL names for vertex and fragment shaders, shader program
	GLuint  vertex;
	GLuint  fragment;
	GLuint  program;

	// initialize shader and program names to zero (OpenGL reserved value)
	MyShader() : vertex(0), fragment(0), program(0)
	{}
};

// load, compile, and link shaders, returning true if successful
bool InitializeShaders(MyShader *shader)
{
	// load shader source from files
	string vertexSource = LoadSource("vertex.glsl");
	string fragmentSource = LoadSource("fragment.glsl");
	if (vertexSource.empty() || fragmentSource.empty()) return false;

	// compile shader source into shader objects
	shader->vertex = CompileShader(GL_VERTEX_SHADER, vertexSource);
	shader->fragment = CompileShader(GL_FRAGMENT_SHADER, fragmentSource);

	// link shader program
	shader->program = LinkProgram(shader->vertex, shader->fragment);

	// check for OpenGL errors and return false if error occurred
	return !CheckGLErrors();
}

// deallocate shader-related objects
void DestroyShaders(MyShader *shader)
{
	// unbind any shader programs and destroy shader objects
	glUseProgram(0);
	glDeleteProgram(shader->program);
	glDeleteShader(shader->vertex);
	glDeleteShader(shader->fragment);
}

// --------------------------------------------------------------------------
// Functions to set up OpenGL buffers for storing geometry data

struct MyGeometry
{
	// OpenGL names for array buffer objects, vertex array object
	GLuint  vertexBuffer;
	GLuint  colourBuffer;
	GLuint  vertexArray;
	GLsizei elementCount;

	// initialize object names to zero (OpenGL reserved value)
	MyGeometry() : vertexBuffer(0), colourBuffer(0), vertexArray(0), elementCount(0)
	{}
};

struct MySphere
{
	glm::vec3 point;
	float radius;
	glm::vec3 colour = glm::vec3(-1.0);
	float reflectivity;
	float phong;

} spheres[maxShapes];

struct MyPlane
{
	glm::vec3 normal;
	glm::vec3 point;
	glm::vec3 colour = glm::vec3(-1.0);
	float reflectivity;
	float phong;

} planes[maxShapes];

struct MyTriangle
{
	glm::vec3 p1;
	glm::vec3 p2;
	glm::vec3 p3;
	glm::vec3 normal;
	glm::vec3 colour = glm::vec3(-1.0);
	float reflectivity;
	float phong;

} triangles[maxShapes];

struct MyLight
{
	glm::vec3 source;
	float intensity;
	float ambient;
	float colour;

} light;

// check ray sphere intersection
glm::vec2 sphereIntersect(MySphere *sphere, glm::vec3 ori, glm::vec3 ray) {
	float radius = sphere->radius;
	glm::vec3 point = sphere->point;
	glm::vec3 ray2 = ori - point;

	// t = (-b +- sqrt(b^2 - 4ac)) / 2a
	float a = dot(ray, ray);
	float b = 2.0 * dot(ray2, ray);
	float c = dot(ray2, ray2) - (radius * radius);

	float t1 = -1.0;
	float t2 = -1.0;

	float det = b*b - 4*a*c;
	if (det > 0.0) {
		float sqRoot = sqrt(det);
		t1 = ((-b) - sqRoot) / (2.0*a);
		t2 = ((-b) + sqRoot) / (2.0*a);
	}

	return glm::vec2(t1, t2);
}

// check ray plane intersection
float planeIntersect(MyPlane *plane, glm::vec3 ori, glm::vec3 ray) {
	glm::vec3 normal = plane->normal;
	float t = -dot(ori - plane->point, normal) / dot(ray, normal);

	if (t > 0.0) return t;
	else return -1.0;
}

// check ray triangle intersection
float triIntersect(MyTriangle *tri, glm::vec3 ori, glm::vec3 ray) {
	glm::vec3 p1 = tri->p1;
	glm::vec3 p2 = tri->p2;
	glm::vec3 p3 = tri->p3;

	glm::vec3 e1 = p2 - p1;
	glm::vec3 e2 = p3 - p1;

	glm::vec3 normal = cross(ray, e2);

	float a = dot(e1, normal);
	if (a > -0.00001 && a < 0.00001) return -1.0;
	float f = 1.0 / a;

	glm::vec3 s = ori - p1;
	float u = f * dot(s, normal);
	if (u < 0.0 || u > 1.0) return -1.0;

	glm::vec3 q = cross(s, e1);
	float v = f * dot(ray, q);
	if (v < 0.0 || u + v > 1.0) return -1.0;

	float t = f * dot(e2, q);
	if (t > 0.0) {
		glm::vec3 point = ori + t * ray;
		tri->normal = cross(p1 - point, p2 - point);
		return t;
	}
	else return -1.0;
}

// get shadow colour
glm::vec3 shadowColour(glm::vec3 point, glm::vec3 viewRay, glm::vec3 colour, glm::vec3 normal, float phong) {
	float t2 = -1.0;
	bool shadows = false;
	glm::vec3 ray2 = glm::normalize(light.source - point);
	float distance = 999.0;
	float d1 = glm::distance(point, light.source);

	// check for intersections
	for (int i = 0; i < maxShapes; i++) {

		// spheres
		if (spheres[i].colour.r != -1.0) {
			glm::vec2 ts = sphereIntersect(&spheres[i], point, ray2);
			if (ts.x > -0.001 && ts.y > -0.001) {
				t2 = ts.x;
				if (t2 > -0.001) {
					float d2 = glm::length(t2 * ray2);
					if (d2 < distance) {
						distance = d2;
					}
				}
			}
		}

		// planes
		if (planes[i].colour.r != -1.0) {
			t2 = planeIntersect(&planes[i], point, ray2);
			if (t2 > 0.0) {
				float d2 = glm::length(t2 * ray2);
				if (d2 < distance && d2 > 0.00001) {
					distance = d2;
				}
			}
		}

		// triangles
		if (triangles[i].colour.r != -1.0) {
			t2 = triIntersect(&triangles[i], point, ray2);
			if (t2 > 0.0) {
				float d2 = glm::length(t2 * ray2);
				if (d2 < distance && d2 > 0.00001) {
					distance = d2;
				}
			}
		}
	}
	if (distance < d1) shadows = true;

	// L = ka*Ia + kd*I*max(0, dot(n,l)) + ks*I*max(0, dot(n,h))^p

	// ka: surface ambient coefficient/"ambient colour"
	// Ia: ambient light intensity
	// kd: diffuse coefficient or surface colour
	// I: intensity of light source
	// n: surface normal
	// l: light vector (light source point - intersection point)
	// ks: specular coefficient
	// h: normalize(v + l)
	// p: phong exponent

	// ambient lighting
	glm::vec3 newColour = colour * light.ambient; // L = ka * Ia

	if (!shadows) {
		// diffuse lighting
		glm::vec3 kd = colour * diffRatio;//colour;
		glm::vec3 norm = glm::normalize(normal);
		glm::vec3 l = glm::normalize(light.source - point);
		newColour += kd * light.intensity * max(0, dot(norm, l));

		// specular lighting
		glm::vec3 h = glm::normalize(l - viewRay);
		float maxTerm = max(0, dot(norm, h));
		float specular = light.colour * light.intensity * pow(maxTerm, phong);
		newColour += glm::vec3(specular);
	}

	return newColour;
}

// get reflection colour
glm::vec3 reflectColour(glm::vec3 ori, glm::vec3 ray, int count) {
	float t = 999.0;
	float reflectivity = 0.0;
	float phong = 0.0;
	glm::vec3 colour = glm::vec3(0.0);
	glm::vec3 normal;

	// check intersections
	for (int i = 0; i < maxShapes; i++) {

		// spheres
		if (spheres[i].colour.r != -1.0) {
			glm::vec2 ts = sphereIntersect(&spheres[i], ori, ray);
			if (ts.x > 0.1 && ts.x < t) {
				t = ts.x;
				colour = spheres[i].colour;
				normal = origin + t * ray - spheres[i].point;
				reflectivity = spheres[i].reflectivity;
				phong = spheres[i].phong;
			}
			else if (ts.y > 0.1 && ts.y < t) {
				t = ts.y;
				colour = spheres[i].colour;
				normal = origin + t * ray - spheres[i].point;
				reflectivity = spheres[i].reflectivity;
				phong = spheres[i].phong;
			}
		}

		// planes
		if (planes[i].colour.r != -1.0) {
			float t1 = planeIntersect(&planes[i], ori, ray);
			if (t1 > 0.1 && t1 < t) {
				t = t1;
				colour = planes[i].colour;
				normal = planes[i].normal;
				reflectivity = planes[i].reflectivity;
				phong = planes[i].phong;
			}
		}

		// triangles
		if (triangles[i].colour.r != -1.0) {
			float t1 = triIntersect(&triangles[i], ori, ray);
			if (t1 > 0.1 && t1 < t) {
				t = t1;
				colour = triangles[i].colour;
				normal = triangles[i].normal;
				reflectivity = triangles[i].reflectivity;
				phong = triangles[i].phong;
			}
		}
	}
	
	t -= 0.1;
	glm::vec3 point = ori + t * ray;
	glm::vec3 newColour = colour;

	if (t > 0.0001 && phong > 0.0) {
		// light and shadow
		newColour = shadowColour(point, ray, colour, normal, phong);

		// reflection
		if (reflectivity > 0.0 && count < 10) {
			glm::vec3 reflectRay = glm::normalize(glm::reflect(ray, glm::normalize(normal)));
			glm::vec3 refColour = reflectColour(point, reflectRay, count + 1);

			newColour.r = (1 - reflectivity) * newColour.r + reflectivity * refColour.r;
			newColour.g = (1 - reflectivity) * newColour.g + reflectivity * refColour.g;
			newColour.b = (1 - reflectivity) * newColour.b + reflectivity * refColour.b;
		}
	}

	return newColour;
}

// get pixel colour
glm::vec3 getColour(int n) {
	glm::vec3 newColour = glm::vec3((3.0 + viewRays[n].x + viewRays[n].y + viewRays[n].z) / 6.0);
	glm::vec3 surfaceNormal;
	float closest = 999.0;
	float reflectivity = 0.0;
	float phong = 0.0;
	float t;
	bool intersect = false;

	// check sphere intersections
	for (int i = 0; i < maxShapes; i++) {
		if (spheres[i].colour.r == -1.0) break;

		glm::vec2 ts = sphereIntersect(&spheres[i], origin, viewRays[n]);
		if (ts.x < 0.0) t = ts.y;
		else t = min(ts.x, ts.y);
		if (t > 0.0 && t < closest) {
			newColour = spheres[i].colour;
			closest = t;
			surfaceNormal = origin + t * viewRays[n] - spheres[i].point;
			reflectivity = spheres[i].reflectivity;
			phong = spheres[i].phong;
			intersect = true;
		}
	}

	//check plane intersections
	for (int i = 0; i < maxShapes; i++) {
		if (planes[i].colour.r == -1.0) break;

		t = planeIntersect(&planes[i], origin, viewRays[n]);
		if (t > 0.0 && t < closest) {
			newColour = planes[i].colour;
			closest = t;
			surfaceNormal = planes[i].normal;
			reflectivity = planes[i].reflectivity;
			phong = planes[i].phong;
			intersect = true;
		}
	}

	// check triangle intersections
	for (int i = 0; i < maxShapes; i++) {
		if (triangles[i].colour.r == -1.0) break;

		t = triIntersect(&triangles[i], origin, viewRays[n]);
		if (t > 0.0 && t < closest) {
			newColour = triangles[i].colour;
			closest = t;
			surfaceNormal = triangles[i].normal;
			reflectivity = triangles[i].reflectivity;
			phong = triangles[i].phong;
			intersect = true;
		}
	}

	if (intersect) {
		t = closest;

		// light and shadows
		glm::vec3 point = origin + t * viewRays[n];
		glm::vec3 L = shadowColour(point, viewRays[n], newColour, surfaceNormal, phong);
		newColour = L;

		// reflections
		if (reflectivity > 0.0) {
			glm::vec3 reflectRay = glm::reflect(viewRays[n], glm::normalize(surfaceNormal));
			glm::vec3 refColour = reflectColour(point, reflectRay, 0);

			newColour.r = (1 - reflectivity) * L.r + reflectivity * refColour.r;
			newColour.g = (1 - reflectivity) * L.g + reflectivity * refColour.g;
			newColour.b = (1 - reflectivity) * L.b + reflectivity * refColour.b;
		}
	}

	// normalize colours
	if (newColour.r < 0.0) newColour.r = 0.0;
	if (newColour.r > 1.0) newColour.r = 1.0;
	if (newColour.g < 0.0) newColour.g = 0.0;
	if (newColour.g > 1.0) newColour.g = 1.0;
	if (newColour.b < 0.0) newColour.b = 0.0;
	if (newColour.b > 1.0) newColour.b = 1.0;

	return newColour;
}

// initialize scene objects
void initializeScene() {

	// clear all shapes
	for (int i = 0; i < maxShapes; i++) {
		spheres[i].colour = glm::vec3(-1.0);
		planes[i].colour = glm::vec3(-1.0);
		triangles[i].colour = glm::vec3(-1.0);
	}

	// scene 1
	if (currScene == 0) {

		// sphere
		spheres[0].radius = 0.825;
		spheres[0].point = glm::vec3(0.9, -1.925, -6.69);
		spheres[0].colour = glm::vec3(0.8, 0.8, 0.8);
		spheres[0].reflectivity = 0.5;
		spheres[0].phong = 1000.0;

		// pyramid
		for (int i = 0; i <= 3; i++) {
			triangles[i].colour = glm::vec3(0.2, 0.8, 0.8);
			triangles[i].reflectivity = 0.3;
			triangles[i].phong = 10.0;
		}

		triangles[0].p1 = glm::vec3(-0.4, -2.75, -9.55);
		triangles[0].p2 = glm::vec3(-0.93, 0.55, -8.51);
		triangles[0].p3 = glm::vec3(0.11, -2.75, -7.98);

		triangles[1].p1 = glm::vec3(0.11, -2.75, -7.98);
		triangles[1].p2 = glm::vec3(-0.93, 0.55, -8.51);
		triangles[1].p3 = glm::vec3(-1.46, -2.75, -7.47);

		triangles[2].p1 = glm::vec3(-1.46, -2.75, -7.47);
		triangles[2].p2 = glm::vec3(-0.93, 0.55, -8.51);
		triangles[2].p3 = glm::vec3(-1.97, -2.75, -9.04);

		triangles[3].p1 = glm::vec3(-1.97, -2.75, -9.04);
		triangles[3].p2 = glm::vec3(-0.93, 0.55, -8.51);
		triangles[3].p3 = glm::vec3(-0.4, -2.75, -9.55);

		// back
		planes[0].normal = glm::vec3(0.0, 0.0, 1.0);
		planes[0].point = glm::vec3(0.0, 0.0, -10.5);
		planes[0].colour = glm::vec3(0.7, 0.7, 0.7);
		planes[0].reflectivity = 0.0;
		planes[0].phong = 2.0;

		// ceiling
		for (int i = 4; i <= 5; i++) {
			triangles[i].colour = glm::vec3(1.0, 1.0, 1.0);
			triangles[i].reflectivity = 0.0;
			triangles[i].phong = 1.0;
		}

		triangles[4].p1 = glm::vec3(2.75, 2.75, -10.5);
		triangles[4].p2 = glm::vec3(2.75, 2.75, -5);
		triangles[4].p3 = glm::vec3(-2.75, 2.75, -5);

		triangles[5].p1 = glm::vec3(-2.75, 2.75, -10.5);
		triangles[5].p2 = glm::vec3(2.75, 2.75, -10.5);
		triangles[5].p3 = glm::vec3(-2.75, 2.75, -5);

		// Floor
		for (int i = 10; i <= 11; i++) {
			triangles[i].colour = glm::vec3(0.6, 0.6, 0.6);
			triangles[i].reflectivity = 0.0;
			triangles[i].phong = 1.0;
		}

		triangles[10].p1 = glm::vec3(2.75, -2.75, -5);
		triangles[10].p2 = glm::vec3(2.75, -2.75, -10.5);
		triangles[10].p3 = glm::vec3(-2.75, -2.75, -10.5);

		triangles[11].p1 = glm::vec3(-2.75, -2.75, -5);
		triangles[11].p2 = glm::vec3(2.75, -2.75, -5);
		triangles[11].p3 = glm::vec3(-2.75, -2.75, -10.5);

		// Green wall on right
		for (int i = 6; i <= 7; i++) {
			triangles[i].colour = glm::vec3(0.2, 0.8, 0.2);
			triangles[i].reflectivity = 0.0;
			triangles[i].phong = 2.0;
		}

		triangles[6].p1 = glm::vec3(2.75, 2.75, -5);
		triangles[6].p2 = glm::vec3(2.75, 2.75, -10.5);
		triangles[6].p3 = glm::vec3(2.75, -2.75, -10.5);

		triangles[7].p1 = glm::vec3(2.75, -2.75, -5);
		triangles[7].p2 = glm::vec3(2.75, 2.75, -5);
		triangles[7].p3 = glm::vec3(2.75, -2.75, -10.5);

		// Red wall on left
		for (int i = 8; i <= 9; i++) {
			triangles[i].colour = glm::vec3(0.8, 0.2, 0.2);
			triangles[i].reflectivity = 0.0;
			triangles[i].phong = 2.0;
		}

		triangles[8].p1 = glm::vec3(-2.75, -2.75, -5);
		triangles[8].p2 = glm::vec3(-2.75, -2.75, -10.5);
		triangles[8].p3 = glm::vec3(-2.75, 2.75, -10.5);

		triangles[9].p1 = glm::vec3(-2.75, 2.75, -5);
		triangles[9].p2 = glm::vec3(-2.75, -2.75, -5);
		triangles[9].p3 = glm::vec3(-2.75, 2.75, -10.5);

		// light
		light.source = glm::vec3(0, 2.5, -7.75);
		light.intensity = 0.65;
		light.ambient = 0.35;
		light.colour = 0.5;
	}

	// scene 2
	else if (currScene == 1) {

		// floor
		planes[0].normal = glm::vec3(0.0, 1.0, 0.0);
		planes[0].point = glm::vec3(0.0, -1.0, 0.0);
		planes[0].colour = glm::vec3(0.8, 0.8, 0.8);
		planes[0].reflectivity = 0.0;
		planes[0].phong = 2.0;

		// back wall
		planes[1].normal = glm::vec3(0.0, 0.0, 1.0);
		planes[1].point = glm::vec3(0.0, 0.0, -12.0);
		planes[1].colour = glm::vec3(0.1, 0.4, 0.4);
		planes[1].reflectivity = 0.0;
		planes[1].phong = 8.0;

		// large yellow sphere
		spheres[0].radius = 0.5;
		spheres[0].point = glm::vec3(1, -0.5, -3.5);
		spheres[0].colour = glm::vec3(0.8, 0.8, 0.2);
		spheres[0].reflectivity = 0.0;
		spheres[0].phong = 30.0;

		// reflective grey sphere
		spheres[1].radius = 0.4;
		spheres[1].point = glm::vec3(0, 1, - 5);
		spheres[1].colour = glm::vec3(0.9, 0.9, 0.9);
		spheres[1].reflectivity = 0.4;
		spheres[1].phong = 200.0;

		// metallic purple sphere
		spheres[2].radius = 0.25;
		spheres[2].point = glm::vec3(-0.8, - 0.75, - 4);
		spheres[2].colour = glm::vec3(0.8, 0.1, 0.8);
		spheres[2].reflectivity = 0.3;
		spheres[2].phong = 120.0;

		// green cone
		for (int i = 0; i <= 11; i++) {
			triangles[i].colour = glm::vec3(0.2, 0.7, 0.2);
			triangles[i].reflectivity = 0.0;
			triangles[i].phong = 80.0;
		}

		triangles[0].p1 = glm::vec3(0, -1, -5.8);
		triangles[0].p2 = glm::vec3(0, 0.6, -5);
		triangles[0].p3 = glm::vec3(0.4, -1, -5.693);

		triangles[1].p1 = glm::vec3(0.4, -1, -5.693);
		triangles[1].p2 = glm::vec3(0, 0.6, -5);
		triangles[1].p3 = glm::vec3(0.6928, -1, -5.4);

		triangles[2].p1 = glm::vec3(0.6928, - 1, - 5.4);
		triangles[2].p2 = glm::vec3(0, 0.6, - 5);
		triangles[2].p3 = glm::vec3(0.8, - 1, - 5);

		triangles[3].p1 = glm::vec3(0.8, - 1, - 5);
		triangles[3].p2 = glm::vec3(0 ,0.6, - 5);
		triangles[3].p3 = glm::vec3(0.6928, - 1, - 4.6);

		triangles[4].p1 = glm::vec3(0.6928, - 1, - 4.6);
		triangles[4].p2 = glm::vec3(0 ,0.6 ,- 5);
		triangles[4].p3 = glm::vec3(0.4 ,- 1, - 4.307);

		triangles[5].p1 = glm::vec3(0.4, - 1, - 4.307);
		triangles[5].p2 = glm::vec3(0, 0.6, - 5);
		triangles[5].p3 = glm::vec3(0, - 1, - 4.2);

		triangles[6].p1 = glm::vec3(0 ,- 1, - 4.2);
		triangles[6].p2 = glm::vec3(0, 0.6, - 5);
		triangles[6].p3 = glm::vec3(-0.4, - 1, - 4.307);

		triangles[7].p1 = glm::vec3(-0.4 ,- 1 ,- 4.307);
		triangles[7].p2 = glm::vec3(0 ,0.6, - 5);
		triangles[7].p3 = glm::vec3(-0.6928 ,- 1, - 4.6);

		triangles[8].p1 = glm::vec3(-0.6928 ,- 1, - 4.6);
		triangles[8].p2 = glm::vec3(0, 0.6 ,- 5);
		triangles[8].p3 = glm::vec3(-0.8, - 1, - 5);

		triangles[9].p1 = glm::vec3(-0.8 ,- 1, - 5);
		triangles[9].p2 = glm::vec3(0 ,0.6 ,- 5);
		triangles[9].p3 = glm::vec3(-0.6928 ,- 1 ,- 5.4);

		triangles[10].p1 = glm::vec3(-0.6928 ,- 1 ,- 5.4);
		triangles[10].p2 = glm::vec3(0 ,0.6 ,- 5);
		triangles[10].p3 = glm::vec3(-0.4, - 1 ,- 5.693);

		triangles[11].p1 = glm::vec3(-0.4 ,- 1 ,- 5.693);
		triangles[11].p2 = glm::vec3(0, 0.6 ,- 5);
		triangles[11].p3 = glm::vec3(0, - 1, - 5.8);

		// Shiny red icosahedron
		for (int i = 12; i <= 31; i++) {
			triangles[i].colour = glm::vec3(1.0, 0.3, 0.3);
			triangles[i].reflectivity = 0.3;
			triangles[i].phong = 100.0;
		}

		triangles[12].p1 = glm::vec3(-2 ,- 1, - 7);
		triangles[12].p2 = glm::vec3(-1.276, - 0.4472, - 6.474);
		triangles[12].p3 = glm::vec3(-2.276, - 0.4472 ,- 6.149);

		triangles[13].p1 = glm::vec3(-1.276 ,- 0.4472, - 6.474);
		triangles[13].p2 = glm::vec3(-2 ,- 1 ,- 7);
		triangles[13].p3 = glm::vec3(-1.276, - 0.4472, - 7.526);

		triangles[14].p1 = glm::vec3(-2 ,- 1 ,- 7);
		triangles[14].p2 = glm::vec3(-2.276 ,- 0.4472, - 6.149);
		triangles[14].p3 = glm::vec3(-2.894, - 0.4472, - 7);

		triangles[15].p1 = glm::vec3(-2 ,- 1, - 7);
		triangles[15].p2 = glm::vec3(-2.894 ,- 0.4472 ,- 7);
		triangles[15].p3 = glm::vec3(-2.276, - 0.4472, - 7.851);

		triangles[16].p1 = glm::vec3(-2 ,- 1 ,- 7);
		triangles[16].p2 = glm::vec3(-2.276, - 0.4472, - 7.851);
		triangles[16].p3 = glm::vec3(-1.276 ,- 0.4472, - 7.526);

		triangles[17].p1 = glm::vec3(-1.276, - 0.4472, - 6.474);
		triangles[17].p2 = glm::vec3(-1.276, - 0.4472, - 7.526);
		triangles[17].p3 = glm::vec3(-1.106 ,0.4472, - 7);

		triangles[18].p1 = glm::vec3(-2.276 ,- 0.4472, - 6.149);
		triangles[18].p2 = glm::vec3(-1.276, - 0.4472, - 6.474);
		triangles[18].p3 = glm::vec3(-1.724 ,0.4472, - 6.149);

		triangles[19].p1 = glm::vec3(-2.894 ,- 0.4472, - 7);
		triangles[19].p2 = glm::vec3(-2.276, - 0.4472, - 6.149);
		triangles[19].p3 = glm::vec3(-2.724 ,0.4472 ,- 6.474);

		triangles[20].p1 = glm::vec3(-2.276 ,- 0.4472, - 7.851);
		triangles[20].p2 = glm::vec3(-2.894, - 0.4472, - 7);
		triangles[20].p3 = glm::vec3(-2.724, 0.4472, - 7.526);

		triangles[21].p1 = glm::vec3(-1.276 ,- 0.4472, - 7.526);
		triangles[21].p2 = glm::vec3(-2.276, - 0.4472, - 7.851);
		triangles[21].p3 = glm::vec3(-1.724, 0.4472, - 7.851);

		triangles[22].p1 = glm::vec3(-1.276 ,- 0.4472, - 6.474);
		triangles[22].p2 = glm::vec3(-1.106 ,0.4472, - 7);
		triangles[22].p3 = glm::vec3(-1.724, 0.4472, - 6.149);

		triangles[23].p1 = glm::vec3(-2.276 ,- 0.4472, - 6.149);
		triangles[23].p2 = glm::vec3(-1.724 ,0.4472, - 6.149);
		triangles[23].p3 = glm::vec3(-2.724, 0.4472, - 6.474);

		triangles[24].p1 = glm::vec3(-2.894 ,- 0.4472, - 7);
		triangles[24].p2 = glm::vec3(-2.724 ,0.4472, - 6.474);
		triangles[24].p3 = glm::vec3(-2.724 ,0.4472, - 7.526);

		triangles[25].p1 = glm::vec3(-2.276 ,- 0.4472, - 7.851);
		triangles[25].p2 = glm::vec3(-2.724, 0.4472 ,- 7.526);
		triangles[25].p3 = glm::vec3(-1.724 ,0.4472, - 7.851);

		triangles[26].p1 = glm::vec3(-1.276, - 0.4472, - 7.526);
		triangles[26].p2 = glm::vec3(-1.724, 0.4472, - 7.851);
		triangles[26].p3 = glm::vec3(-1.106 ,0.4472, - 7);

		triangles[27].p1 = glm::vec3(-1.724 ,0.4472 ,- 6.149);
		triangles[27].p2 = glm::vec3(-1.106 ,0.4472 ,- 7);
		triangles[27].p3 = glm::vec3(-2, 1 ,- 7);

		triangles[28].p1 = glm::vec3(-2.724, 0.4472, - 6.474);
		triangles[28].p2 = glm::vec3(-1.724, 0.4472 ,- 6.149);
		triangles[28].p3 = glm::vec3(-2 ,1 ,- 7);

		triangles[29].p1 = glm::vec3(-2.724, 0.4472, - 7.526);
		triangles[29].p2 = glm::vec3(-2.724 ,0.4472 ,- 6.474);
		triangles[29].p3 = glm::vec3(-2, 1, - 7);

		triangles[30].p1 = glm::vec3(-1.724, 0.4472, - 7.851);
		triangles[30].p2 = glm::vec3(-2.724, 0.4472 ,- 7.526);
		triangles[30].p3 = glm::vec3(-2, 1, - 7);

		triangles[31].p1 = glm::vec3(-1.106, 0.4472, - 7);
		triangles[31].p2 = glm::vec3(-1.724, 0.4472, - 7.851);
		triangles[31].p3 = glm::vec3(-2, 1 ,- 7);

		// light
		light.source = glm::vec3(4, 6, -1);
		light.intensity = 0.9;
		light.ambient = 0.3;
		light.colour = 0.5;
	}

	// scene 3
	else if (currScene == 2) {

		// back walls
		for (int i = 0; i <= 3; i++) {
			planes[i].colour = glm::vec3(0.0, 0.2, 0.3);
			planes[i].reflectivity = 0.3;
			planes[i].phong = 3.0;
		}

		planes[0].normal = glm::vec3(-0.005, -0.005, 1.0);
		planes[0].point = glm::vec3(0.0, 0.0, -20.0);

		planes[1].normal = glm::vec3(-0.005, 0.005, 1.0);
		planes[1].point = glm::vec3(0.0, 0.0, -20.0);

		planes[2].normal = glm::vec3(0.005, -0.005, 1.0);
		planes[2].point = glm::vec3(0.0, 0.0, -20.0);

		planes[3].normal = glm::vec3(-0.005, -0.005, 1.0);
		planes[3].point = glm::vec3(0.0, 0.0, -20.0);

		// spiral of spheres and triangles
		for (int i = 0; i < maxShapes; i++) {
			float div = 3.0 * maxShapes / 32.0;
			float t = float(i) / div;
			float phong = float(i + 1) / float(maxShapes + 1);

			// spheres
			spheres[i].point = glm::vec3(0.8 * t * sin(t), 
										0.8 * t * cos(t), 
										-8.0 - t);
			spheres[i].radius = t / 8.0;
			spheres[i].colour = glm::vec3(sin(t)*sin(t), cos(t)*cos(t), sin(t)*cos(t));
			spheres[i].reflectivity = 0.8 * float(i) / float(maxShapes);
			spheres[i].phong = 2000.0 * phong;

			// triangles
			t += 2.0 / div;
			float diff = (1.0 / div) / 3.0;
			triangles[i].p1 = glm::vec3(-0.8 * t * sin(t), 
										-0.8 * t * cos(t),
										-8.0 - t);
			triangles[i].p2 = glm::vec3(-0.8 * (t + diff) * sin((t + diff))+sin(t), 
										-0.8 * (t + diff) * cos((t + diff))+cos(t), 
										-8.0 - (t + diff));
			triangles[i].p3 = glm::vec3(-0.8 * (t + 2*diff) * sin((t + 2 * diff)), 
										-0.8 * (t + 2 * diff) * cos((t + 2 * diff)), 
										-8.0 - (t + 2 * diff));
			triangles[i].colour = glm::vec3(sin(t)*sin(t), cos(t)*cos(t), sin(t)*cos(t));
			triangles[i].reflectivity = 0.0;
			triangles[i].phong = 100.0 * phong;
		}

		// light
		light.source = glm::vec3(0, 0, -4);
		light.intensity = 1.0;
		light.ambient = 0.7;
		light.colour = 0.5;
	}
}

// create buffers and fill with geometry data, returning true if successful
bool InitializeGeometry(MyGeometry *geometry)
{
	// initialize objects for current scene
	initializeScene();

	GLfloat z = -sqrt(wWidth*wWidth + wHeight*wHeight) / (4.0 * tan(fov / 2.0));

	// initialize vertices/view rays and set colours
	for (int i = 0; i < wWidth*wHeight; i++) {
		GLfloat x = (float(i % wWidth) - float(wWidth) / 2.0);
		GLfloat y = (float(i / wWidth) - float(wHeight) / 2.0);

		vertices[i][0] = 2.0 * (x) / float(wWidth);
		vertices[i][1] = 2.0 * (y) / float(wHeight);

		viewRays[i] = glm::normalize(glm::vec3(x, y, z));

		glm::vec3 newColour = getColour(i);
		colours[i][0] = newColour.r;
		colours[i][1] = newColour.g;
		colours[i][2] = newColour.b;
	}
	
	geometry->elementCount = wWidth * wHeight;

	// these vertex attribute indices correspond to those specified for the
	// input variables in the vertex shader
	const GLuint VERTEX_INDEX = 0;
	const GLuint COLOUR_INDEX = 1;

	// create an array buffer object for storing our vertices
	glGenBuffers(1, &geometry->vertexBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, geometry->vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	// create another one for storing our colours
	glGenBuffers(1, &geometry->colourBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, geometry->colourBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(colours), colours, GL_STATIC_DRAW);

	// create a vertex array object encapsulating all our vertex attributes
	glGenVertexArrays(1, &geometry->vertexArray);
	glBindVertexArray(geometry->vertexArray);

	// associate the position array with the vertex array object
	glBindBuffer(GL_ARRAY_BUFFER, geometry->vertexBuffer);
	glVertexAttribPointer(VERTEX_INDEX, 2, GL_FLOAT, GL_FALSE, 0, 0); //2
	glEnableVertexAttribArray(VERTEX_INDEX);

	// assocaite the colour array with the vertex array object
	glBindBuffer(GL_ARRAY_BUFFER, geometry->colourBuffer);
	glVertexAttribPointer(COLOUR_INDEX, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(COLOUR_INDEX);

	// unbind our buffers, resetting to default state
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	// check for OpenGL errors and return false if error occurred
	return !CheckGLErrors();
}

// deallocate geometry-related objects
void DestroyGeometry(MyGeometry *geometry)
{
	// unbind and destroy our vertex array object and associated buffers
	glBindVertexArray(0);
	glDeleteVertexArrays(1, &geometry->vertexArray);
	glDeleteBuffers(1, &geometry->vertexBuffer);
	glDeleteBuffers(1, &geometry->colourBuffer);
}

// --------------------------------------------------------------------------
// Rendering function that draws our scene to the frame buffer

void RenderScene(MyGeometry *geometry, MyShader *shader)
{
	// bind our shader program and the vertex array object containing our
	// scene geometry, then tell OpenGL to draw our geometry
	glUseProgram(shader->program);
	glBindVertexArray(geometry->vertexArray);
	glDrawArrays(GL_POINTS, 0, geometry->elementCount);

	// reset state to default (no shader or geometry bound)
	glBindVertexArray(0);
	glUseProgram(0);

	// check for an report any OpenGL errors
	CheckGLErrors();
}

// --------------------------------------------------------------------------
// GLFW callback functions

// reports GLFW errors
void ErrorCallback(int error, const char* description)
{
	cout << "GLFW ERROR " << error << ":" << endl;
	cout << description << endl;
}

// handles keyboard input events
void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

	// set scene 1
	else if (key == GLFW_KEY_1 && action == GLFW_RELEASE) {
		if (currScene != 0) {
			currScene = 0;
			initializeScene();
			initialized = false;
		}
	}

	// set scene 2
	else if (key == GLFW_KEY_2 && action == GLFW_RELEASE) {
		if (currScene != 1) {
			currScene = 1;
			initializeScene();
			initialized = false;
		}
	}

	// set scene 3
	else if (key == GLFW_KEY_3 && action == GLFW_RELEASE) {
		if (currScene != 2) {
			currScene = 2;
			initializeScene();
			initialized = false;
		}
	}

	// save image
	else if (key == GLFW_KEY_S && action == GLFW_RELEASE) {
		ImageWriter image_writer(wWidth, wHeight);
		for (int x = 0; x < wWidth; x++) {
			for (int y = 0; y < wHeight; y++) {
				image_writer.set(x, wHeight - y - 1,
					int(float(255) * colours[wWidth*y + x][0]),
					int(float(255) * colours[wWidth*y + x][1]),
					int(float(255) * colours[wWidth*y + x][2]));
			}
		}
		string fname = "scene" + std::to_string(currScene + 1) + ".png";
		char *fchar = const_cast<char*>(fname.c_str());
		image_writer.save(fchar);
	}
}

// ==========================================================================
// PROGRAM ENTRY POINT

int main(int argc, char *argv[])
{
	// initialize the GLFW windowing system
	if (!glfwInit()) {
		cout << "ERROR: GLFW failed to initialize, TERMINATING" << endl;
		return -1;
	}
	glfwSetErrorCallback(ErrorCallback);

	// attempt to create a window with an OpenGL 4.1 core profile context
	GLFWwindow *window = 0;
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	if (antialiasing) glfwWindowHint(GLFW_SAMPLES, 4); // antialiasing

	window = glfwCreateWindow(wWidth, wHeight, "CPSC 453 A4", 0, 0);
	if (!window) {
		cout << "Program failed to create GLFW window, TERMINATING" << endl;
		glfwTerminate();
		return -1;
	}

	// set keyboard callback function and make our context current (active)
	glfwSetKeyCallback(window, KeyCallback);
	glfwMakeContextCurrent(window);

	//Intialize GLAD
#ifndef LAB_LINUX
	if (!gladLoadGL())
	{
		cout << "GLAD init failed" << endl;
		return -1;
	}
#endif

	// query and print out information about our OpenGL environment
	QueryGLVersion();

	// call function to load and compile shader programs
	MyShader shader;
	if (!InitializeShaders(&shader)) {
		cout << "Program could not initialize shaders, TERMINATING" << endl;
		return -1;
	}

	// call function to create and fill buffers with geometry data
	MyGeometry geometry;
	

	// run an event-triggered main loop
	while (!glfwWindowShouldClose(window))
	{
		if (!initialized) {
			DestroyGeometry(&geometry);
			if (!InitializeGeometry(&geometry))
				cout << "Program failed to intialize geometry!" << endl;
			initialized = true;
		}

		RenderScene(&geometry, &shader);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	// clean up allocated resources before exit
	DestroyGeometry(&geometry);
	DestroyShaders(&shader);
	glfwDestroyWindow(window);
	glfwTerminate();

	cout << "Goodbye!" << endl;
	return 0;
}

// ==========================================================================
// SUPPORT FUNCTION DEFINITIONS

// --------------------------------------------------------------------------
// OpenGL utility functions

void QueryGLVersion()
{
	// query opengl version and renderer information
	string version = reinterpret_cast<const char *>(glGetString(GL_VERSION));
	string glslver = reinterpret_cast<const char *>(glGetString(GL_SHADING_LANGUAGE_VERSION));
	string renderer = reinterpret_cast<const char *>(glGetString(GL_RENDERER));

	cout << "OpenGL [ " << version << " ] "
		<< "with GLSL [ " << glslver << " ] "
		<< "on renderer [ " << renderer << " ]" << endl;
}

bool CheckGLErrors()
{
	bool error = false;
	for (GLenum flag = glGetError(); flag != GL_NO_ERROR; flag = glGetError())
	{
		cout << "OpenGL ERROR:  ";
		switch (flag) {
		case GL_INVALID_ENUM:
			cout << "GL_INVALID_ENUM" << endl; break;
		case GL_INVALID_VALUE:
			cout << "GL_INVALID_VALUE" << endl; break;
		case GL_INVALID_OPERATION:
			cout << "GL_INVALID_OPERATION" << endl; break;
		case GL_INVALID_FRAMEBUFFER_OPERATION:
			cout << "GL_INVALID_FRAMEBUFFER_OPERATION" << endl; break;
		case GL_OUT_OF_MEMORY:
			cout << "GL_OUT_OF_MEMORY" << endl; break;
		default:
			cout << "[unknown error code]" << endl;
		}
		error = true;
	}
	return error;
}

// --------------------------------------------------------------------------
// OpenGL shader support functions

// reads a text file with the given name into a string
string LoadSource(const string &filename)
{
	string source;

	ifstream input(filename.c_str());
	if (input) {
		copy(istreambuf_iterator<char>(input),
			istreambuf_iterator<char>(),
			back_inserter(source));
		input.close();
	}
	else {
		cout << "ERROR: Could not load shader source from file "
			<< filename << endl;
	}

	return source;
}

// creates and returns a shader object compiled from the given source
GLuint CompileShader(GLenum shaderType, const string &source)
{
	// allocate shader object name
	GLuint shaderObject = glCreateShader(shaderType);

	// try compiling the source as a shader of the given type
	const GLchar *source_ptr = source.c_str();
	glShaderSource(shaderObject, 1, &source_ptr, 0);
	glCompileShader(shaderObject);

	// retrieve compile status
	GLint status;
	glGetShaderiv(shaderObject, GL_COMPILE_STATUS, &status);
	if (status == GL_FALSE)
	{
		GLint length;
		glGetShaderiv(shaderObject, GL_INFO_LOG_LENGTH, &length);
		string info(length, ' ');
		glGetShaderInfoLog(shaderObject, info.length(), &length, &info[0]);
		cout << "ERROR compiling shader:" << endl << endl;
		cout << source << endl;
		cout << info << endl;
	}

	return shaderObject;
}

// creates and returns a program object linked from vertex and fragment shaders
GLuint LinkProgram(GLuint vertexShader, GLuint fragmentShader)
{
	// allocate program object name
	GLuint programObject = glCreateProgram();

	// attach provided shader objects to this program
	if (vertexShader)   glAttachShader(programObject, vertexShader);
	if (fragmentShader) glAttachShader(programObject, fragmentShader);

	// try linking the program with given attachments
	glLinkProgram(programObject);

	// retrieve link status
	GLint status;
	glGetProgramiv(programObject, GL_LINK_STATUS, &status);
	if (status == GL_FALSE)
	{
		GLint length;
		glGetProgramiv(programObject, GL_INFO_LOG_LENGTH, &length);
		string info(length, ' ');
		glGetProgramInfoLog(programObject, info.length(), &length, &info[0]);
		cout << "ERROR linking shader program:" << endl;
		cout << info << endl;
	}

	return programObject;
}
