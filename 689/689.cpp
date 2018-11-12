
#include "Field.h"
#include <iostream>
#include "Camera.h"
#include "Renderer.h"
#include "Implicit.h"
#include "Transforms.h"
#include "Grids.h"
#include "CSG.h"
#include "CPT.h"
#include <random>
#include <list>
//#include "OBJ.h"
#include "OBJ2.h"
#include <glm/glm.hpp>

using namespace vmr;

int main()
{
	ColorField gray = ColorField(Color(0.7, 0.7, 0.7, 1));

	Light light1 = Light(Vector(-3, -3, -3), Color(1, 1, 1, 1));
	Light light2 = Light(Vector(2, 0, 1), Color(0.5, 0.5, 0.5, 1));
	Light light3 = Light(Vector(0, -0.5, -2), Color(0.5, 0.5, 0.5, 1));

	Camera c = Camera();
	c.setEyeViewUp(Vector(0,0,5), Vector(0, 0, -1), Vector(0, 1, 0));
	c.setNearPlane(0.5);
	c.setFarPlane(15);

	char objFile[] = "bunnyLow.obj";
	char file1[] = "bunnyGrid025.bin";
	char file2[] = "bunnyGridGrad1.bin";
	char file3[] = "bunnyGridF1.bin";
	char fileName[] = "bun.exr";
	FieldFloatPtr lightTmp[] = { NULL };

	std::vector < glm::vec3 > vertices;
	std::vector < glm::vec2 > uvs;
	std::vector < glm::vec3 > normals;

	Polygon poly;
	poly.loadOBJ(objFile, vertices, uvs, normals);
	LevelSetFieldFaces bunFaces = LevelSetFieldFaces(poly.getFaces());
	stampDensityGrid(Vector(-1, -1, -1), Vector(2, 2, 2), 0.025, &bunFaces, file3);


	//loadOBJ(objFile, vertices, uvs, normals);
	//LevelSetGradField bunGrad = LevelSetGradField(&vertices,&normals);
	//stampVectorGrid(Vector(-1, -1, -1), Vector(2, 2, 2), 0.1, &bunGrad, file2);
	//float * gradArr = loadVectorGrid(file2);
	//griddedVectorField bunGradG = griddedVectorField(gradArr, Vector(-1, -1, -1), Vector(2, 2, 2), 0.1);
	//LevelSetField bun = LevelSetField(&vertices, &normals);
	//stampDensityGrid(Vector(-1, -1, -1), Vector(2, 2, 2), 0.1, &bun, file1);
	
	float * arr = loadGrid(file3);
	griddedFloatField bunF = griddedFloatField(arr,Vector(-1, -1, -1), Vector(2, 2, 2), 0.025);
	FloatScale flip = FloatScale(&bunF, Vector(1, -1, 1));
	Sphere topArt = Sphere(Vector(0, -1.73, 0),1);
	Plane botArt = Plane(Vector(0, 0.75, 0), Vector(0, 1, 0));
	Sphere earArt = Sphere(Vector(-1.1, -0.5, -0.5), 0.5);
	Cutout bun1 = Cutout(&flip, &topArt);
	Cutout bun2 = Cutout(&bun1, &botArt);
	Cutout bun3 = Cutout(&bun2, &earArt);
	//CPTLevel bunCPT = CPTLevel(&bun3, &bunGradG);
	NoiseScalar bunNoise = NoiseScalar(3.1, 2.5, 0.6, 0.1, 0.5, 2);
	//Cloud bunC = Cloud(&bun3, &bunNoise, &bunCPT);

	Sphere earsOnly = Sphere(Vector(-0.25, -1.1, 0), 0.75);
	GradVectorField vec = GradVectorField(&topArt);
	FSPNVectorField vec2 = FSPNVectorField(&topArt,&bunNoise,0.1);
	CutoutVec vecEars = CutoutVec(&vec2, &earsOnly);
	FloatAdvect ad = FloatAdvect(&bun3, &vecEars, 2);
	
	//Renderer r = Renderer(&c, 1920, 1080, 0.01, &bun1, &gray, lightTmp); 
	Renderer r = Renderer(&c, 1920, 1080, 0.01, &bun3, &gray, lightTmp); 
	//Renderer r = Renderer(&c2, 640, 360, 0.01, &ns14, &gray, lightTmp); 
	r.render(fileName);


}
