#pragma once
#ifndef VMR_GRIDS_H
#define VMR_GRIDS_H

#include "Field.h"
#include "Vector.h"
#include <iostream>
#include <fstream>
#include "OBJ.h"

namespace vmr
{

	void stampDensityGrid(Vector llc, Vector n,float step, FieldFloatPtr floatField, char * fileName)
	{
		std::ofstream myfile;
		myfile.open(fileName);
		myfile << llc.X() << " " << llc.Y() << " " << llc.Z() << "\n";
		myfile << n.X() << " " << n.Y() << " " << n.Z() << "\n";
		myfile << step << "\n";

		int width = 50;
		float * arr = NULL;
		int numPoints = width*width*width;
		arr = new float[numPoints];
		
		int count = 0;
		for (int k = 0; k < width; k++)
		{
			for (int j = 0; j < width; j++)
			{
#pragma omp parallel for 
				for (int i = 0; i < width; i++)
				{
					Vector currentPoint = llc + Vector(i*step, j*step, k*step);
					//calculate density at point
					float d = floatField->eval(currentPoint);
					//assign to grid
					//myfile << d << " ";
					int index = (i + width * (j + width * k));
					arr[index] = d;
					count++;
					std::cout << "\n" << count;
				}
			}
		}
		for (int i = 0; i < numPoints; i++)
		{
			myfile << arr[i] << " ";
		}
		myfile.close();
	}

	void stampWispGrid(Vector llc, Vector n, float step, Noise_t t1, Noise_t t2, Wisp w, char * fileName)
	{
		std::ofstream myfile;
		myfile.open(fileName, std::ios::binary);
		myfile << llc.X() << " " << llc.Y() << " " << llc.Z() << "\n";
		myfile << n.X() << " " << n.Y() << " " << n.Z() << "\n";
		myfile << step << "\n";

		int width = floor(n[0] / step);
		float * arr = NULL;
		int numPoints = width*width*width;
		arr = new float[numPoints];

		std::mt19937 gen1(123);
		std::uniform_real_distribution<double> rand(-1.0, 1.0);

		int count = 0;
		#pragma omp parallel for 
		for (int i = 0; i < w.numDots; i++)
		{
			//generate dot
			
			//place in random position in box
			Vector r0 = Vector(rand(gen1), rand(gen1), rand(gen1));
			//move to unit sphere
			Vector r1 = r0.unitvector();
			//move off unit sphere randomly
			FractalSum<PerlinNoiseGustavson> fs;
			fs.setParameters(t1);
			float fsum1 = fs.eval(r0);
			//std::cout << "\n fsum1: " << fsum1;
			float q = pow(abs(fsum1), w.clump);
			Vector r2 = r1 * q;
			//put in object space
			Vector P2 = w.pos + (r2 * w.pscale);
			//add displacement
			fs.setParameters(t2);
			float fsum2 = fs.eval(P2);
			//std::cout << "\n fsum2: " << fsum1;
			float fsum3 = fs.eval(P2 + w.offset);
			float fsum4 = fs.eval(P2 - w.offset);
			Vector D = Vector(fsum2, fsum3, fsum4);
			Vector P3 = P2 + (w.offScale*D);
			//evaluate position in grid  || reverse trilinear interpolation
			Vector xyz = Vector(floor(((P3.X() - llc.X()) / step)), floor(((P3.Y() - llc.Y()) / step)), floor(((P3.Z() - llc.Z()) / step)));
			int ii = xyz.X();
			int j = xyz.Y();
			int k = xyz.Z();
			int t = floor(P3.X() / step);
			float wi = (P3.X() - t * step) / step;
			t = floor(P3.Y() / step);
			float wj = (P3.Y() - t * step) / step;
			t = floor(P3.Z() / step);
			float wk = clip((P3.Z() - t * step) / step, 0.00001, 1);

			int c000 = (ii + width * (j + width * k));
			int c100 = (ii + 1) + width * (j + width * k);
			int c001 = (ii)+width * ((j + 1) + width * k);
			int c010 = (ii)+width * (j + width * (k + 1));
			int c101 = (ii + 1) + width * ((j + 1) + width * (k));
			int c110 = (ii + 1) + width * (j + width * (k + 1));
			int c011 = (ii)+width * ((j + 1) + width * (k + 1));
			int c111 = (ii + 1) + width * ((j + 1) + width * (k + 1));

			arr[c000] += w.density * (1.0 - wi) * (1.0 - wj) * (1.0 - wk);
			arr[c100] += w.density * wi * (1.0 - wj) * (1.0 - wk);
			arr[c001] += w.density * (1.0 - wi) * wj * (1.0 - wk);
			arr[c010] += w.density * (1.0 - wi) * (1.0 - wj) * (wk);
			arr[c101] += w.density * (wi) * (wj) * (1.0 - wk);
			arr[c110] += w.density * (wi) * (1.0 - wj) * (wk);
			arr[c011] += w.density * (1.0 - wi) * (wj) * (wk);
			arr[c111] += w.density * wi * wj * wk;

			//std::cout << "\n dot num: " << i;
		}
		std::cout << "\nfinished dots";

		for (int i = 0; i < numPoints; i++)
		{
			myfile << arr[i] << " ";
		}
		myfile.close();
	}

	void stampWispGridArr(Vector llc, Vector n, float step, Noise_t t1, Noise_t t2, Wisp w, float * &array)
	{
		int width = floor(n[0] / step);
		float * arr = NULL;
		int numPoints = width * width*width;
		arr = new float[numPoints];

		std::mt19937 gen1(123);
		std::uniform_real_distribution<double> rand(-1.0, 1.0);

		int count = 0;
		#pragma omp parallel for 
		for (int i = 0; i < w.numDots; i++)
		{
			//generate dot

			//place in random position in box
			Vector r0 = Vector(rand(gen1), rand(gen1), rand(gen1));
			//move to unit sphere
			Vector r1 = r0.unitvector();
			//move off unit sphere randomly
			FractalSum<PerlinNoiseGustavson> fs;
			fs.setParameters(t1);
			float fsum1 = fs.eval(r0);
			float q = pow(abs(fsum1), w.clump);
			Vector r2 = r1 * q;
			//put in object space
			Vector P2 = w.pos + (r2 * w.pscale);
			//add displacement
			fs.setParameters(t2);
			float fsum2 = fs.eval(P2);
			float fsum3 = fs.eval(P2 + w.offset);
			float fsum4 = fs.eval(P2 - w.offset);
			Vector D = Vector(fsum2, fsum3, fsum4);
			Vector P3 = P2 + (w.offScale*D);
			//evaluate position in grid  || reverse trilinear interpolation
			if (isIn(P3, n, (llc + n / 2)))
			{
				float ii = (P3[0] - llc[0]) / step;
				float j = (P3[1] - llc[1]) / step;
				float k = (P3[2] - llc[2]) / step;

				float wi = ii - (int)(ii);
				float wj = j - (int)(j);
				float wk = k - (int)(k);

				ii = (int)(ii);
				j = (int)(j);
				k = (int)(k);

				if (ii < width - 1 && j < width - 1 && k < width - 1)
				{
					int c000 = (ii + width * (j + width * k));
					int c100 = (ii + 1) + width * (j + width * k);
					int c001 = (ii)+width * ((j + 1) + width * k);
					int c010 = (ii)+width * (j + width * (k + 1));
					int c101 = (ii + 1) + width * ((j + 1) + width * (k));
					int c110 = (ii + 1) + width * (j + width * (k + 1));
					int c011 = (ii)+width * ((j + 1) + width * (k + 1));
					int c111 = (ii + 1) + width * ((j + 1) + width * (k + 1));

					arr[c000] += w.density * (1.0 - wi) * (1.0 - wj) * (1.0 - wk);
					arr[c100] += w.density * wi * (1.0 - wj) * (1.0 - wk);
					arr[c001] += w.density * (1.0 - wi) * wj * (1.0 - wk);
					arr[c010] += w.density * (1.0 - wi) * (1.0 - wj) * (wk);
					arr[c101] += w.density * (wi) * (wj) * (1.0 - wk);
					arr[c110] += w.density * (wi) * (1.0 - wj) * (wk);
					arr[c011] += w.density * (1.0 - wi) * (wj) * (wk);
					arr[c111] += w.density * wi * wj * wk;
				}
			}
			//std::cout << "\n dot num: " << i;
		}
		array = arr;
	}

	void stampLightGrid(Vector llc, Vector n, float step, FieldFloatPtr floatField, Light light, char * fileName)
	{
		std::ofstream myfile;
		myfile.open(fileName);
		myfile << llc.X() << " " << llc.Y() << " " << llc.Z() << "\n";
		myfile << n.X() << " " << n.Y() << " " << n.Z() << "\n";
		myfile << step << "\n";
		float * arr = NULL;
		int numPoints = (n.X() / step)*(n.Y() / step)*(n.Z() / step);
		arr = new float[numPoints];
		int width = ceil(n.Z() / step);
		int count = 0;
		for (int k = 0; k < width; k++)
		{
			for (int j = 0; j < width; j++)
			{
				#pragma omp parallel for 
				for (int i = 0; i < width; i++)
				{
					Vector currentPoint = llc + Vector(i*step, j*step, k*step);
					//calculate density at point
					float tL = 0;
					Vector lightRay = (light.location() - currentPoint).unitvector();
					Vector currentPosL = currentPoint;
					float dsm = 0;
					float ds = 0.01;
					int iter = (light.location() - currentPoint).magnitude() / ds;
					for (int k = 0; k < iter; k++)
					{
						currentPosL = currentPosL + (lightRay * ds);
						dsm += mask(floatField, currentPosL) * ds;
					}
					tL = std::exp(-4 * dsm);
					//assign to grid
					//myfile << tL << " ";
					int index = (i + width * (j + width * k));
					arr[index] = tL;
					count++;
					std::cout << "\n" << count;
				}
			}
		}
		for (int i = 0; i < numPoints; i++)
		{
			myfile << arr[i] << " ";
		}
		myfile.close();
	}

	float * loadGrid(char * fileName)
	{
		std::ifstream file(fileName, std::ios::binary);
		float * arr = NULL;
		if (file.is_open())
		{
			float llcX,llcY,llcZ,nX,nY,nZ,step;
			file >> llcX;
			file >> llcY;
			file >> llcZ;
			file >> nX;
			file >> nY;
			file >> nZ;
			file >> step;
			int numPoints = (nX / step)*(nY / step)*(nZ / step);
			arr = new float[numPoints];
			for (int i = 0; i < numPoints; i++)
			{
				file >> arr[i];
			}
		}
		return arr;
	}

}

#endif VMR_GRIDS_H