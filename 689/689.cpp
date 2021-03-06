
#include "Field.h"
#include <iostream>
#include "Camera.h"
#include "Renderer.h"
#include "Implicit.h"
#include "Transforms.h"
#include "Grids.h"
#include "CSG.h"
#include <random>
#include <list>
#include "OBJ.h"
#include <glm/glm.hpp>

using namespace vmr;

int main()
{
	ColorField gray = ColorField(Color(0.7, 0.7, 0.7, 1));

	Light light1 = Light(Vector(-1, -2, 2), Color(1, 1, 1, 1));
	Light light2 = Light(Vector(2, 0, 1), Color(0.5, 0.5, 0.5, 1));
	Light light3 = Light(Vector(0, -0.5, -2), Color(0.5, 0.5, 0.5, 1));

	Camera c = Camera();
	c.setEyeViewUp(Vector(0, 0, 5), Vector(0, 0, -1), Vector(0, 1, 0));
	c.setNearPlane(3);
	c.setFarPlane(15);

	char fileName1[] = "testWispGrid.bin";
	

	FieldFloatPtr lightTmp[] = { NULL };

	Wisp w = {
		Vector(0,0,0), //pos
		1, //pscale
		1, //density
		1, //clump
		Vector(0.5,0.5,0.5), //offset
		0.5, //offScale
		5000000 //numDots
	};

	Noise_t t1;
	t1.amplitude = 1;
	t1.fjump = 2.0;
	t1.frequency = 2.0;
	t1.octaves = 2;
	t1.pscale = 1;
	Noise_t t2;
	t2.amplitude = 1;
	t2.fjump = 2.5;
	t2.frequency = 3.5;
	t2.octaves = 3;
	t2.pscale = 1;

	char fileName[] = "wisp.exr";

	//std::cout << "\nstamping wisp grid: ";
	//stampWispGrid(Vector(-2.5, -2.5, -2.5), Vector(5, 5, 5), 0.01, t1, t2, w, fileName1);
	//std::cout << "\nloading wisp grid: ";
	//float * arr = loadGrid(fileName1);
	//float * arr;
	//stampWispGridArr(Vector(-1.5, -1.5, -1.5), Vector(3, 3, 3), 0.01, t1, t2, w, arr);
	//griddedFloatField wispF = griddedFloatField(arr, Vector(-1.5, -1.5, -1.5), Vector(3, 3, 3), 0.01);
	//std::cout << "\n rendering . . . ";
	//Renderer r2 = Renderer(&c, 960, 540, 0.005, &wispF, &gray, lightTmp);
	//r2.render(fileName);

	float clump = 1;
	float o = 1;
	float f = 0.5;
	float fJ = 0.5;
	int frame = 1;
	for (int count = 1; count < (500 - frame); count++)
	{
		std::cout << "\n now rendering frame: " << frame + (count-1);
		std::cout << "\n with o=" << o << " f=" << f << " fJ=" << fJ << " clump=" << clump;
						
		float * arr;
		stampWispGridArr(Vector(-1.5, -1.5, -1.5), Vector(3, 3, 3), 0.006, t1, t2, w, arr);
		griddedFloatField wispF = griddedFloatField(arr, Vector(-1.5, -1.5, -1.5), Vector(3, 3, 3), 0.006);

		char integer_string[4];
		char name[50] = "wisp_";
		const char *ext = ".exr";
		sprintf_s(integer_string, "%d", frame + (count - 1));
		strcat_s(name, integer_string);
		strcat_s(name, ext);
		Renderer r2 = Renderer(&c, 960, 540, 0.001, &wispF, &gray, lightTmp);
		r2.render(name);

		if (f + 1 < 5) 
			f++;
		else
		{
			f = 0.5;
			if (fJ + 1 < 4) 
				f++;
			else
			{
				fJ = 0.5;
				if (o + 0.5 < 5)
					o += 0.5;
				else
				{
					o = 1;
					if (clump + 0.5 < 4) 
						clump += 0.5;
					else
					{
						break;
					}
				}
			}
		}
			
	}

	//NoiseStamp ns = NoiseStamp(Vector(0, 0, 0), 1, 3, 1.1, 1, 1, 3, 2);
	//char name[50] = "noise.exr";
	//Renderer r2 = Renderer(&c, 960, 540, 0.01, &ns, &gray, lightTmp);
	//r2.render(name);

	//int count = 1;
	//for (float o = 1; o < 8; o += 1) // start at 1
	//{
	//	for (float fa = 1; fa < 3; fa += 0.33)// start at 1
	//	{
	//		for (float fJump = 0.5; fJump < 4; fJump += 1) // start at 0.5
	//		{
	//			for (float f = 0.5; f < 5; f += 1)// start at 0.5
	//			{
	//				std::cout << "\n now rendering frame: " << count;
	//				std::cout << "\n with o=" << o << " f=" << f << " fJ=" << fJump << " fa=" << fa;
	//				
	//				NoiseStamp ns = NoiseStamp(Vector(0, 0, 0), 1, f, fJump, 0.6, 1, o, fa);

	//				char integer_string[4];
	//				char name[50] = "noiseStamp_";
	//				const char *ext = ".exr";
	//				sprintf_s(integer_string, "%d", count);
	//				strcat_s(name, integer_string);
	//				strcat_s(name, ext);

	//				Renderer r2 = Renderer(&c, 960, 540, 0.001, &ns, &gray, lightTmp);
	//				//Renderer r2 = Renderer(&c, 640, 360, 0.001, &ns, &gray, lightTmp);
	//				r2.render(name);

	//				count++;
	//			}
	//		}
	//	}
	//}

	//PyroSphere ps = PyroSphere(Vector(0, 0, 0), 0.5, 2, 1.5, 0.6, 1, 0.3, 3);
	//char name[50] = "pyro_test.exr";
	//Renderer r2 = Renderer(&c, 960, 540, 0.01, &ps, &gray, lightTmp);
	//Renderer r2 = Renderer(&c, 640, 360, 0.01, &ps, &gray, lightTmp);
	//r2.render(name);

	//int count = 1;
	//for (float o = 1; o < 4; o += 1) // start at 1
	//{
	//	for (float g = 0.01; g < 2; g += 0.225)// start at .01
	//	{
	//		for (float f = 0.5; f < 5; f += 1) // start at .5
	//		{
	//			for (float fJump = 0.5; fJump < 4; fJump += 1)// start at .5
	//			{
	//				std::cout << "\n now rendering frame: " << count;
	//				std::cout << "\n with o=" << o << " f=" << f << " fJ=" << fJump << " g=" << g;
	//				
	//				PyroSphere ps = PyroSphere(Vector(0, 0, 0), 0.5, f, fJump, 0.6, 1, g, o);

	//				char integer_string[4];
	//				char name[50] = "pyroSphere_";
	//				const char *ext = ".exr";
	//				sprintf_s(integer_string, "%d", count);
	//				strcat_s(name, integer_string);
	//				strcat_s(name, ext);

	//				Renderer r2 = Renderer(&c, 960, 540, 0.001, &ps, &gray, lightTmp);
	//				r2.render(name);

	//				count++;
	//			}
	//		}
	//	}
	//}
}
