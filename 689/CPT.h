#pragma once
#ifndef CPT_H
#define CPT_H

#include "Vector.h"
#include "Field.h"
#include "PerlinNoise1.h"
#include "Noise.h"

namespace vmr
{
	class TerrainNoise : public Field<float>
	{

	public:

		TerrainNoise(Noise_t t2, float aUp, float aDown, float gammaUp, float gammaDown)
		{
			t1 = t2;
			au = aUp;
			ad = aDown;
			gu = gammaUp;
			gd = gammaDown;
		}

		const float eval(const Vector& P) const
		{
			FractalSum<PerlinNoiseGustavson> fs;
			fs.setParameters(t1);
			float n = fs.eval(P);
			if (n >= 0)
				return au * pow(n, gu);
			else
				return -(ad*pow(abs(n), gd));
		}

	private:

		Noise_t t1;
		float au;
		float ad;
		float gu;
		float gd;
	};

	class NoiseScalar : public Field<float>
	{

	public:

		NoiseScalar(float freq, float fJump, float roughness, float A, float gamm, float oct)
		{
			t.gamma = gamm;
			t.frequency = freq;
			t.fjump = fJump;
			t.octaves = oct;
			t.amplitude = A;
			t.roughness = roughness;
		}

		const float eval(const Vector& P) const
		{

			FractalSum<PerlinNoiseGustavson> fs;
			fs.setParameters(t);
			float n = fs.eval(P);
			return t.amplitude * pow(abs(n), t.gamma);
		}

	private:
		Noise_t t;
	};

	class CPTImplicit : public Field<Vector>
	{

	public:

		CPTImplicit(Field<float> * f1)
		{
			f = f1;
		}

		const Vector eval(const Vector& P) const
		{
			return (P - (f->eval(P)*f->grad(P)));
		}

	private:
		FieldFloatPtr f;
	};

	class CPTLevel : public Field<Vector>
	{

	public:

		CPTLevel(Field<float> * f1,Field<Vector> * g1)
		{
			f = f1;
			g = g1;
		}

		const Vector eval(const Vector& P) const
		{
			return (P - (f->eval(P)*g->eval(P).unitvector()));
		}

	private:
		FieldFloatPtr f;
		FieldVectorPtr g;
	};

	class NPTImplicit : public Field<Vector>
	{

	public:

		NPTImplicit(Field<float> * f1)
		{
			f = f1;
		}

		const Vector eval(const Vector& P) const
		{
			Vector x = P - f->eval(P) * (f->grad(P) / pow(f->grad(P).magnitude(), 2));
			while (f->eval(x) < -0.0001 || f->eval(x) > 0.0001)
			{
				x = x - f->eval(x) * (f->grad(x) / pow(f->grad(x).magnitude(), 2));
			}
			return x;
		}

	private:
		FieldFloatPtr f;
	};

	class Terrain : public Field<float>
	{

	public:

		Terrain(Field<float> * f1, FieldFloatPtr f2)
		{
			f = f1;
			g = f2;
		}

		const float eval(const Vector& P) const
		{
			if (f->eval(P) + g->eval(P) >= 0)
				return 1;
			else
				return 0;
			//return f->eval(P) + g->eval(P);
		}

	private:
		FieldFloatPtr f;
		FieldFloatPtr g;
	};

	class Warp : public Field<float>
	{

	public:

		Warp(Field<float> * f1, Field<Vector> * v1)
		{
			f = f1;
			v = v1;
		}

		const float eval(const Vector& P) const
		{
			return f->eval(v->eval(P));
		}

	private:
		FieldFloatPtr f;
		FieldVectorPtr v;
	};

	class Cloud : public Field<float>
	{

	public:

		Cloud(Field<float> * f1, Field<float> * noise, Field<Vector> * CPT)
		{
			f = f1;
			n = noise;
			v = CPT;
		}

		const float eval(const Vector& P) const
		{
			return f->eval(P) + (n->eval(P)*v->eval(P)).magnitude();
		}

	private:
		FieldFloatPtr f;
		FieldFloatPtr n;
		FieldVectorPtr v;
	};



}

#endif // !CPT_H
