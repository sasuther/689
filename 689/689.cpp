
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
#include "OBJ.h"
#include <glm/glm.hpp>

using namespace vmr;

int main()
{
	ColorField gray = ColorField(Color(0.7, 0.7, 0.7, 1));

	Light light1 = Light(Vector(-3, -3, -3), Color(1, 1, 1, 1));
	Light light2 = Light(Vector(2, 0, 1), Color(0.5, 0.5, 0.5, 1));
	Light light3 = Light(Vector(0, -0.5, -2), Color(0.5, 0.5, 0.5, 1));

	Camera c = Camera();
	c.setEyeViewUp(Vector(0.2,-1, 2), Vector(0, 0.3, -1), Vector(0, 1, 0));
	c.setNearPlane(0.5);
	c.setFarPlane(15);

	Camera c2 = Camera();
	c2.setEyeViewUp(Vector(2, -1, 2), Vector(-0.5, 0.3, -0.5), Vector(0, 1, 0));
	c2.setNearPlane(0.5);
	c2.setFarPlane(15);

	Camera c3 = Camera();
	c3.setEyeViewUp(Vector(-2, -1, 0), Vector(1, 0.3, 0), Vector(0, 1, 0));
	c3.setNearPlane(0.5);
	c3.setFarPlane(15);

	Camera c4 = Camera();
	c4.setEyeViewUp(Vector(0, 0, 7), Vector(0, 0, -1), Vector(0, 1, 0));
	c4.setNearPlane(0.5);
	c4.setFarPlane(15);

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
	t2.amplitude = 0.5;
	t2.fjump = 2.1;
	t2.frequency = 2.5;
	t2.octaves = 3;
	t2.pscale = 1;

	Plane terPlane = Plane(Vector(0, 0, 0), Vector(0, 1, 0));
	CPTImplicit terCPT = CPTImplicit(&terPlane);
	TerrainNoise terNoise = TerrainNoise(t1, 1, 0, 1, 0.2);
	Warp terWarp = Warp(&terNoise, &terCPT);
	Terrain terrain = Terrain(&terPlane, &terWarp);
	Box terBox = Box(Vector(0, -0.5, 0), 7, 6);
	Intersection terInt = Intersection(&terrain, &terBox);

	char fileName2[] = "terrain.bin";
	char fileName3[] = "output.csv";
	char fileName4[] = "noiseGrid.bin";
	char fileName5[] = "teapot_grid.txt";
	//stampDensityGrid(Vector(-3.5, -4, -3.5), Vector(7,7,7), 0.05, &terInt, fileName2);
	//std::cout << "\n Load density grid";
	//float * terArray = loadGrid(fileName2);
	//griddedFloatField terGrid = griddedFloatField(terArray, Vector(-3.5, -4, -3.5), Vector(7, 7, 7), 0.05);
	//std::cout << "\n Stamp light grid";
	//stampLightGrid(Vector(-3.5, -4, -3.5), Vector(7, 7, 7), 0.1, &terGrid, light1, fileName3);
	//std::cout << "\n Load light grid";
	//float * terLightArr = loadGrid(fileName3);
	//griddedFloatField terLights = griddedFloatField(terLightArr, Vector(-3.5, -4, -3.5), Vector(7, 7, 7), 0.1);

	//FieldFloatPtr lightTmp[] = { &terLights };

	/*NPTImplicit terNPT = NPTImplicit(&terPlane);
	Warp terWarp2 = Warp(&terNoise, &terNPT);
	Terrain terrain2 = Terrain(&terPlane, &terWarp2);
	Box terBox2 = Box(Vector(0, -0.5, 0), 7, 6);
	Intersection terInt2 = Intersection(&terrain2, &terBox2);*/

	Torus bridge = Torus(Vector(-0.3, 0.2, 0.45), 0.45, 0.1, Vector(0.2, 0, -1));
	CPTImplicit bridgeCPT = CPTImplicit(&bridge);
	NoiseScalar bridgeNoise = NoiseScalar(3.5,1.5,0.6,0.015,0.5,4);
	Warp bridgeWarp = Warp(&bridgeNoise, &bridgeCPT);
	Terrain bridgeEnd = Terrain(&bridge, &bridgeWarp);

	PyroSphere volcano = PyroSphere(Vector(-0.625, -0.55, .225), 0.2,5.2,1.3,1,0.01,0.77,4);

	Union end = Union(&terInt, &bridgeEnd);
	//Union end1 = Union(&end, &volcano);
	Cutout end1 = Cutout(&end, &volcano);

	PyroSphere cave = PyroSphere(Vector(0.55, -0.15, -0.2), 0.15,4.3,2.1,0.9,0.01,0.6,3);
	//Sphere cave = Sphere(Vector(0.55, -0.1, -0.2), 0.15);
	Cutout end2 = Cutout(&end1, &cave);

	//float * pArr = loadParticles(fileName3);
	//std::cout << "\n stamping . . .";
	//stampNoiseGrid(Vector(-3, -3, -3), Vector(6, 6, 6), 0.05, partArr, 26, fileName4);
	//float * noiseArr = loadGrid(fileName4);
	//griddedFloatField smoke = griddedFloatField(noiseArr, Vector(-3, -3, -3), Vector(6, 6, 6), 0.05);
	NoiseStamp * nsA = new NoiseStamp[26];
	for (int i = 0; i < 26; i++)
	{
		//Vector pos = Vector(pArr[i]-0.6, pArr[i + 1]-0.3, pArr[i + 2]);
		Vector pos = Vector(-0.65, -0.35 - (i*0.025), 0.23 - (i*0.03));
		nsA[i] = NoiseStamp(pos, 0.05 + (i * 0.01), 4.1, 3.5, 0.75, 0.2, 3, 2 + (i * 0.05));
	}
	//nsA[0] = NoiseStamp(Vector(0, 0, 0), 0.4, 4.1, 3.5, 0.75, 0.2, 4, 3);
	//nsA[1] = NoiseStamp(Vector(1, 0, 0), 0.6, 3.5, 4.0, 0.75, 0.2, 4, 4);
	Union ns1 = Union(&nsA[0], &nsA[1]);
	Union ns2 = Union(&ns1, &nsA[2]);
	Union ns3 = Union(&ns2, &nsA[3]);
	Union ns4 = Union(&ns3, &nsA[4]);
	Union ns5 = Union(&ns4, &nsA[5]);
	Union ns6 = Union(&ns5, &nsA[6]);
	Union ns7 = Union(&ns6, &nsA[7]);
	Union ns8 = Union(&ns7, &nsA[8]);
	Union ns9 = Union(&ns8, &nsA[9]);
	Union ns10 = Union(&ns9, &nsA[10]);
	Union ns11 = Union(&ns10, &nsA[11]);
	Union ns12 = Union(&ns11, &nsA[12]);
	Union ns13 = Union(&ns12, &nsA[13]);
	Union ns14 = Union(&ns13, &nsA[14]);
	Union ns15 = Union(&ns14, &nsA[15]);
	Union ns16 = Union(&ns15, &nsA[16]);
	Union ns17 = Union(&ns16, &nsA[17]);
	Union ns18 = Union(&ns17, &nsA[18]);
	Union ns19 = Union(&ns18, &nsA[19]);
	Union ns20 = Union(&ns19, &nsA[20]);
	Union ns21 = Union(&ns20, &nsA[21]);
	Union ns22 = Union(&ns21, &nsA[22]);
	Union ns23 = Union(&ns22, &nsA[23]);
	Union ns24 = Union(&ns23, &nsA[24]);
	Union ns25 = Union(&ns24, &nsA[25]);

	Union end3 = Union(&end2, &ns25);

	float * arr = loadGrid(fileName5);
	griddedFloatField mon = griddedFloatField(arr, Vector(-1.5, -1.5, -1.5), Vector(3, 3, 3), 0.1);
	FloatScale monS = FloatScale(&mon, Vector(1, -1, 1));
	Plane monP1 = Plane(Vector(0, -0.65, 0), Vector(0, 1, 0));
	Plane monP2 = Plane(Vector(1.65, 0, 0), Vector(-1, 0, 0));
	//Intersection int2 = Intersection(&int1, &monP2);
	FloatRotation monRot = FloatRotation(&monS, Vector(0, 1, 0), 90);
	FloatScale monS2 = FloatScale(&monRot, Vector(0.05, 0.05, 0.05));
	FloatTranslate monT = FloatTranslate(&monS2, Vector(-1.55, -.75, 0.1));
	Intersection int1 = Intersection(&monT, &monP1);
	
	Union end4 = Union(&end3, &monT);


	char fileName[] = "terrain2.exr";

	Renderer r = Renderer(&c3, 192, 108, 0.01, &end3, &gray, lightTmp); 
	//Renderer r = Renderer(&c2, 640, 360, 0.01, &ns14, &gray, lightTmp); 
	r.render(fileName);


}
