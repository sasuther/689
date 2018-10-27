#pragma once

#ifndef VMR_IMPLICIT_H
#define VMR_IMPLICIT_H

#include "Field.h"
#include "Vector.h"
#include "FSPN.h"
#include "Noise.h"
#include "PerlinNoise1.h"
#include <random>

namespace vmr
{

	class Sphere : public Field<float>
	{

	public:
		
		Sphere(Vector center, float radius)
		{
			cen = center;
			rad = radius;
		}

		const float eval(const Vector& P) const
		{
			return rad - (P - cen).magnitude();
		}
		const Vector grad(const Vector& P) const
		{
			return (-(P - cen) / (P - cen).magnitude());
		}

	private:

		Vector cen;
		float rad;
		
	};

	class PyroSphere : public Field<float>
	{

	public:

		PyroSphere(Vector center, float radius, float freq, float fJump, float roughness,float A,float gamm, float oct)
		{
			t.gamma = gamm;
			t.frequency = freq;
			t.fjump = fJump;
			t.octaves = oct;
			t.amplitude = A;
			t.roughness = roughness;
			t.P = center;
			t.pscale = radius;
		}

		const float eval(const Vector& P) const
		{

			FractalSum<PerlinNoiseGustavson> fs;
			fs.setParameters(t);
			float sum = fs.eval((P-t.P).unitvector()*t.pscale);
			float n = t.amplitude * pow(abs(sum), t.gamma);
			if ((P - t.P).magnitude() < t.pscale + n)
				return 1;
			else
				return 0;

		}

	private:

		Noise_t t;
	
	};

	class NoiseStamp : public Field<float>
	{

	public:

		NoiseStamp(Vector pos, float pscale, float freq, float fJump, float roughness, float A, int oct,float fade)
		{
			fad = fade;
			t.frequency = freq;
			t.fjump = fJump;
			t.octaves = oct;
			t.amplitude = A;
			t.roughness = roughness;
			t.P = pos;
			t.pscale = pscale;
		}

		const float eval(const Vector& P) const
		{
			if ((P - t.P).magnitude() <= t.pscale)
			{

				FractalSum<PerlinNoiseGustavson> fs;
				fs.setParameters(t);
				float nv = fs.eval(P);
				float ff = pow((1 - ((P - t.P).magnitude() / t.pscale)), fad);
				nv *= ff;

				if (nv > 0)
					return nv;
				else
					return 0;
			}
			return 0;
		}

	private:
		float fad;
		Noise_t t;
	};

	class Torus : public Field<float>
	{

	public:

		Torus(Vector center, float radiusMajor, float radiusMinor, Vector axis)
		{
			cen = center;
			radMajor = radiusMajor;
			radMinor = radiusMinor;
			n = axis;
		}

		const float eval(const Vector& P) const
		{
			Vector x = (P - cen);
			Vector xPerp = x - (x*n.unitvector())*n.unitvector();
			return (4 * (radMajor*radMajor)*(xPerp.magnitude()*xPerp.magnitude())) - pow((double)((x.magnitude()*x.magnitude()) + (radMajor*radMajor) - (radMinor*radMinor)), 2.0);
		}
		const Vector grad(const Vector& P) const
		{
			return Vector();
		}

	private:

		Vector cen;
		float radMajor;
		float radMinor;
		Vector n;

	};

	class Cone : public Field<float>
	{

	public:

		Cone(Vector center, float radius, float height, Vector axis)
		{
			cen = center;
			rad = radius;
			h = height;
			n = axis;
			theta = 2 * asin(rad / (sqrt(h*h + rad * rad)));
		}

		const float eval(const Vector& P) const
		{
			if (P == cen)
			{
				return 0;
			}
			else if ((P - cen)*n.unitvector() > h)
			{
				return (h - (P - cen)*n.unitvector());
			}
			else if ((P - cen)*n.unitvector() < 0)
			{
				return (P - cen)*n.unitvector();
			}
			else
			{
				return theta - acos(((P - cen)*n.unitvector()) / (P - cen).magnitude());
			}
			
		}
		const Vector grad(const Vector& P) const
		{
			return Vector();
		}

	private:

		Vector cen;
		float rad;
		float h;
		Vector n;
		float theta;

	};

	class Plane : public Field<float>
	{

	public:

		Plane(Vector center, Vector axis)
		{
			cen = center;
			n = axis.unitvector();
		}

		const float eval(const Vector& P) const
		{
			return (P - cen)*n;
		}
		const Vector grad(const Vector& P) const
		{
			return Vector();
		}

	private:

		Vector cen;
		Vector n;

	};

	class Box : public Field<float> 
	{

	public:

		Box(Vector center, float radius, float pNorm)
		{
			cen = center;
			rad = radius;
			pN = pNorm;
		}

		const float eval(const Vector& P) const
		{
			Vector x = (P - cen);
			return (rad - pow(x[0], pN) - pow(x[1], pN) - pow(x[2], pN));
		}
		const Vector grad(const Vector& P) const
		{
			return Vector();
		}

	private:

		Vector cen;
		float rad;
		float pN;

	};
	
	class Icosahedron : public Field<float> {

		const double PI = 3.14159265;
		const double T = 1.6180339;

	public:
		Icosahedron(Vector center) {
			cen = center;
		}

		const float eval(const Vector&P)const {
			float xMag = (P - cen).magnitude();
			float x = (P - cen).X();
			float y = (P - cen).Y();
			float z = (P - cen).Z();

			if (xMag > (PI*1.8))
				return -1.8*PI;
			else
				return cos(x + T * y) + cos(x - T * y)
				+ cos(y + T * z) + cos(y - T * z)
				+ cos(z + T * x) + cos(z - T * x)
				- 2;
		}

		const Vector grad(const Vector&P)const {
			return Vector();
		}

	private:
		Vector cen;
	};

	class SteinerPatch : public Field<float> 
	{

	public:

		SteinerPatch(Vector center, Vector axis)
		{
			cen = center;
			n = axis.unitvector();
		}

		const float eval(const Vector& P) const
		{
			Vector x = (P - cen);
			return -((x.X()*x.X()*x.Y()*x.Y()) + (x.X()*x.X()*x.Z()*x.Z()) + (x.Y()*x.Y()*x.Z()*x.Z()) - (x.X()*x.Y()*x.Z()));
		}
		const Vector grad(const Vector& P) const
		{
			return Vector();
		}

	private:

		Vector cen;
		Vector n;

	};

	class Cylinder : public Field<float>
	{

	public:

		Cylinder(Vector center, float radius, Vector axis)
		{
			cen = center;
			rad = radius;
			n = axis.unitvector();
		}

		const float eval(const Vector& P) const
		{
			Vector x = (P - cen);
			return (rad - (x - (x*n)*n).magnitude());
		}
		const Vector grad(const Vector& P) const
		{
			return Vector();
		}

	private:

		Vector cen;
		Vector n;
		float rad;
	};

	class Ellipse : public Field<float>
	{

	public:

		Ellipse(Vector center, Vector axis, float radiusMajor, float radiusMinor)
		{
			cen = center;
			n = axis.unitvector();
			radMajor = radiusMajor;
			radMinor = radiusMinor;
		}

		const float eval(const Vector& P) const
		{
			float z = (P - cen)*n;
			Vector xPerp = (P - cen) - z * n;
			return 1 - ((z*z) / (radMajor*radMajor)) - ((xPerp.magnitude()*xPerp.magnitude()) / (radMinor*radMinor));
		}
		const Vector grad(const Vector& P) const
		{
			return Vector();
		}

	private:

		Vector cen;
		Vector n;
		float radMajor;
		float radMinor;
	};

	

	class Object
	{
		public:

			Object(FieldFloatPtr shape, Color color, float density, float kappa)
			{
				s = shape;
				c = color;
				d = density;
				k = kappa;
			}

			float density() 
			{ 
				return d; 
			}
			float kappa() { return k; }
			Color color() { return c; }
			FieldFloatPtr shape()
			{
				return s;
			}

		private:

			FieldFloatPtr s;
			Color c;
			float d;
			float k;
	};

}

#endif // VMR_IMPLICIT_H
