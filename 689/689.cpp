
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
#include <iostream>
#include <iomanip>
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

	Camera c2 = Camera();
	c2.setEyeViewUp(Vector(5, 0, 0), Vector(-1, 0, 0), Vector(0, 1, 0));
	c2.setNearPlane(0.5);
	c2.setFarPlane(15);

	char objFile[] = "bunnyLow.obj";
	char file1[] = "bunnyGrid025.bin";
	char file2[] = "bunnyGridGrad1.bin";
	char file3[] = "bunnyGridF1.bin";
	char fileName[] = "bun.exr";
	FieldFloatPtr lightTmp[] = { NULL };

	std::vector < glm::vec3 > vertices;
	std::vector < glm::vec2 > uvs;
	std::vector < glm::vec3 > normals;

	//Polygon poly;
	//poly.loadOBJ(objFile, vertices, uvs, normals);
	//LevelSetField bunFaces = LevelSetField(poly.getFaces());
	//stampDensityGrid(Vector(-1, -1, -1), Vector(2, 2, 2), 0.05, &bunFaces, file3);
	
	float * arr = loadGrid(file3);
	griddedFloatField bunF = griddedFloatField(arr,Vector(-1, -1, -1), Vector(2, 2, 2), 0.05);
	FloatScale flip = FloatScale(&bunF, Vector(1, -1, 1));
	Sphere topArt = Sphere(Vector(0, -1.73, 0),1);
	Plane botArt = Plane(Vector(0, 0.75, 0), Vector(0, 1, 0));
	Plane pl = Plane(Vector(0, 0, -0.6), Vector(0, 0, -1));
	Sphere earArt = Sphere(Vector(-1.1, -0.5, -0.5), 0.5);
	Box box = Box(Vector(0.8, -1.24, -1), 0.001, 20);
	FloatScale bS = FloatScale(&box, Vector(0.5, 0.5, 0.5));
	Cutout bun1 = Cutout(&flip, &topArt);
	Cutout bun2 = Cutout(&bun1, &botArt);
	Cutout bun3 = Cutout(&bun2, &earArt);
	Cutout bun4 = Cutout(&bun3, &bS);
	Cutout bun5 = Cutout(&bun4, &pl);
	Sphere centerSphere = Sphere(Vector(0, 0, 0), 1);
	//CPTLevel bunCPT = CPTLevel(&bun3, &bunGradG);
	//NPTImplicit bunCPT = NPTImplicit(&bun5);
	CPTImplicit bunCPT = CPTImplicit(&bun5);
	NoiseScalar bunNoise = NoiseScalar(3.5, 2.5, 0.6, 0.2, 0.5, 3);
	Cloud bunC = Cloud(&bun5, &bunNoise, &bunCPT);
	//FSPNVectorField vec1 = FSPNVectorField(&bunC, &bunNoise, 0.1);
	GradVectorField vec = GradVectorField(&centerSphere);
	FloatAdvectItr ad = FloatAdvectItr(&bunC, &vec, 0.003,1);

	/*Sphere earsOnly = Sphere(Vector(-0.25, -1.1, 0), 0.75);
	GradVectorField vec = GradVectorField(&topArt);
	FSPNVectorField vec2 = FSPNVectorField(&topArt,&bunNoise,0.1);
	CutoutVec vecEars = CutoutVec(&vec2, &earsOnly);
	FloatAdvect ad = FloatAdvect(&bun5, &vecEars, 2);*/
	
	//Renderer r = Renderer(&c, 1920, 1080, 0.01, &bun1, &gray, lightTmp); 
	//Renderer r = Renderer(&c, 640, 360, 0.01, &ad, &gray, lightTmp); 
	//Renderer r = Renderer(&c2, 640, 360, 0.01, &ns14, &gray, lightTmp); 
	//r.render(fileName);

	std::stringstream file_2;
	file_2 << "bunCloud_" << std::setw(4) << std::setfill('0') << 1 << ".exr";


	//render frame with no change
	Renderer r = Renderer(&c, 1920, 1080, 0.01, &bun5, &gray, lightTmp);
	r.render((char *)file_2.str().c_str());
	
	//render increasing fluff
	for (int i = 2; i < 62; i++)
	{
		std::stringstream file;
		file << "bunCloud_" << std::setw(4) << std::setfill('0') << i << ".exr";

		bunNoise = NoiseScalar(3.5, 2.5, 0.6, (0.01 + (0.033*(i-2))), 0.5, 3);
		bunC = Cloud(&bunF, &bunNoise, &bunCPT);
		Renderer r = Renderer(&c, 1920, 1080, 0.01, &bunC, &gray, lightTmp);
		r.render((char *)file.str().c_str());
	}

	//render increasing advection

	for (int i = 62; i < 121; i++)
	{
		std::stringstream file;
		file << "bunCloud2_" << std::setw(4) << std::setfill('0') << i << ".exr";
		ad = FloatAdvectItr(&bunC, &vec, 0.003,(i-62)+1);
		
		Renderer r = Renderer(&c, 1920, 1080, 0.01, &ad, &gray, lightTmp);
		r.render((char *)file.str().c_str());
		std::cout << "\nrendered frame: " << i;
	}
}
