#pragma once

#ifndef RENDERER_H
#define RENDERER_H

#include "Camera.h"
#include "Vector.h"
#include "Implicit.h"
#include "Classes.h"
#include <cmath>
#include <random>
#include <iterator>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <thread>
#include "OBJ.h"

#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"

namespace vmr
{
	float mask(FieldFloatPtr f, Vector &x);
	float clamp(FieldFloatPtr f, float a, float b, Vector &x);

	void exportImage(int width, int height, float arr[], char * file);

	float lightCheck(Vector currentPoint, Vector lightLoc, FieldFloatPtr floatField);

	bool isInBounds(Vector pos, Vector boxDim, Vector boxCent);

	class Renderer
	{
		public:

			//Renderer(Camera camera, int imageWidth, int imageHeight, float Sdelta, Object arr[],int numObjects,char * file)
			Renderer(Camera * camera, int imageWidth, int imageHeight, float Sdelta, FieldFloatPtr floatField, FieldColorPtr colField, FieldFloatPtr * lightArr)
			{
				c = camera;
				width1 = imageWidth;
				height1 = imageHeight;
				deltaS1 = Sdelta;
				floatF = floatField;
				colorF = colField;
				lightA = lightArr;
				//lights = light;
				//numLights = numLight;
			}

			void render(char * fileName)
			{
				float* arr = new float[width1*height1 * 4];
				int ofs = 0;
				std::mt19937 generator(123);
				std::uniform_real_distribution<double> dis(0.0, 1.0);
				int count = 0;
				int start = 0;

				for (int j = 0; j < (int)height1; j++)
				{
#pragma omp parallel for 
					for (int i = 0; i < (int)width1; i++)
					{
							int start = (i * 4) + (j*width1 * 4);

							Color L;
							float T;
							float deltaT;

							//for each pixel calculate normal of ray
							float tmp = tan(2 * 3.14159*(c->fov() / 2) / 360);
							float tmp2 = (float)i / (width1 - 1.0f);
							float u = (-1.0f + (dis(generator) / (width1 - 1)) + (2.0f * tmp2)) * tmp;

							tmp = (tan(2 * 3.14159*(c->fov() / 2) / 360) / c->aspectRatio());
							tmp2 = (float)j / (height1 - 1.0f);
							float v = (-1.0f + (dis(generator) / (width1 - 1)) + (2.0f*tmp2)) * tmp;

							Vector q = (u*(c->up() ^ c->view())) + (v*c->up());
							Vector rayNormal = (q + c->view()) / (q + c->view()).magnitude();
							//initialize position on ray
							Vector currentPos = c->eye() + (rayNormal*c->nearPlane());
							L = Color(0, 0, 0, 0);
							T = 1;
							int iterations = (c->farPlane() - c->nearPlane()) / deltaS1;

							deltaT = 0;

							//MARCH!
							for (int t = 0; t < iterations; t++)
							{
								//if in bounding box, compute
								currentPos += rayNormal.unitvector() * deltaS1;
								Color c = Color(0, 0, 0, 0);
								float d = 0;

								if(isInBounds(currentPos, Vector(1.5,1.5,1.5), Vector(0, 0, 0)))
								{
									d = mask(floatF, currentPos);
									//d = clamp(floatF,0,1,currentPos);
									//d = floatF->eval(currentPos);
									
									//loop through lights
									for (int m = 0; m < 1; m++)
									{
										float tL = 0;
										Color lightColor = Color(1, 1, 1, 1);
										Color lightColor2 = Color(0.5, 0.5, 0.5, 1);
										//if (m == 1)
										//	lightColor = Color(0.5, 0.5, 0.5, 1);

										if (d != 0)
										{
											//c = Color(1, 1, 1, 1);
											//c += (lightColor * Color(1,1,1,1) * lightA[m]->eval(currentPos));
											c += (lightColor * Color(1, 1, 1, 1) * lightCheck(currentPos, Vector(-1, -1, 1), floatF));
											//c += (lightColor2 * Color(1, 1, 1, 1) * lightCheck(currentPos, Vector(2, -2, 0), floatF));
										}
					
									}
									deltaT = std::exp(-10 * deltaS1*d);
									L += c * (1 - deltaT)*T;
									T *= deltaT;
								}

							}
								
							arr[start] = L.X();
							arr[start + 1] = L.Y();
							arr[start + 2] = L.Z();
							arr[start + 3] = (1.0f - T);
					}
					//std::cout << "\nline: " << j;
				}
				
				exportImage(width1, height1, arr, fileName);
			}

		private:
			Camera * c;
			int width1;
			int height1;
			float deltaS1;
			FieldFloatPtr floatF;
			FieldColorPtr colorF;
			float deltaT;
			Color L;
			float T;
			//int numLights;
			FieldFloatPtr * lightA;

	};

	float mask(FieldFloatPtr f, Vector &x)
	{
		if (f->eval(x) > 0)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	float clamp(FieldFloatPtr f,float a, float b, Vector &x)
	{ 
		float tmp = f->eval(x);
		if (tmp >= a && tmp <=b)
		{
			return tmp;
		}
		else if(tmp < a)
		{
			return a;
		}
		else
		{
			return b;
		}
	}

	bool isInBounds(Vector pos, Vector boxDim, Vector boxCent)
	{
		if (pos.X() >= (boxCent.X() - boxDim.X() / 2) && pos.X() <= (boxCent.X() + boxDim.X() / 2))
		{
			if (pos.Y() >= (boxCent.Y() - boxDim.Y() / 2) && pos.Y() <= (boxCent.Y() + boxDim.Y() / 2))
			{
				if (pos.Z() >= (boxCent.Z() - boxDim.Z() / 2) && pos.Z() <= (boxCent.Z() + boxDim.Z() / 2))
				{
					return true;
				}
			}
		}
		return false;
	}

	float lightCheck(Vector currentPoint, Vector lightLoc, FieldFloatPtr floatField)
	{
		float tL = 0;
		Vector lightRay = (lightLoc - currentPoint).unitvector();
		Vector currentPosL = currentPoint;
		float dsm = 0;
		float ds = 0.05;
		int iter = (lightLoc - currentPoint).magnitude() / ds;
		for (int k = 0; k < iter; k++)
		{
			currentPosL = currentPosL + (lightRay * ds);
			dsm += mask(floatField, currentPosL) * ds;
		}
		tL = std::exp(-10 * dsm);
		return tL;
	}

	void exportImage(int width, int height,float arr[], char * file)
	{
		EXRHeader header;
		InitEXRHeader(&header);

		EXRImage image;
		InitEXRImage(&image);

		image.num_channels = 4;

		std::vector<float> images[4];
		images[0].resize(width * height);
		images[1].resize(width * height);
		images[2].resize(width * height);
		images[3].resize(width * height);

		int ofs = 0;
		for (int i = 0; i < width * height; i++) {
			images[0][i] = arr[ofs++];
			images[1][i] = arr[ofs++];
			images[2][i] = arr[ofs++];
			images[3][i] = arr[ofs++];
			
		}

		float* image_ptr[4];
		image_ptr[0] = &(images[3].at(0)); // A
		image_ptr[1] = &(images[2].at(0)); // B
		image_ptr[2] = &(images[1].at(0)); // G
		image_ptr[3] = &(images[0].at(0)); // R

		image.images = (unsigned char**)image_ptr;
		image.width = width;
		image.height = height;

		header.num_channels = 4;
		header.channels = (EXRChannelInfo *)malloc(sizeof(EXRChannelInfo) * header.num_channels);
		// Must be BGR(A) order, since most of EXR viewers expect this channel order.
		strncpy_s(header.channels[0].name, "B", 255); header.channels[0].name[strlen("B")] = '\0';
		strncpy_s(header.channels[1].name, "G", 255); header.channels[1].name[strlen("G")] = '\0';
		strncpy_s(header.channels[2].name, "R", 255); header.channels[2].name[strlen("R")] = '\0';
		strncpy_s(header.channels[3].name, "A", 255); header.channels[3].name[strlen("A")] = '\0';

		header.pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
		header.requested_pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
		for (int i = 0; i < header.num_channels; i++) {
			header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
			header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of output image to be stored in .EXR
		}

		const char* err;
		int ret = SaveEXRImageToFile(&image, &header, file, &err);
		if (ret != TINYEXR_SUCCESS) {
			fprintf(stderr, "Save EXR err: %s\n", err);
		}
		else
		{
			//std::cout << "file made";
		}

		free(header.channels);
		free(header.pixel_types);
		free(header.requested_pixel_types);
	}

}

#endif // !RENDERER_H
