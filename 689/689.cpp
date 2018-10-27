
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

	char fileName[] = "wisp.exr";


}
