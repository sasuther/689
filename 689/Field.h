
#pragma once 

#ifndef FIELD_H
#define FIELD_H

#include "Vector.h"
#include "Matrix.h"
#include "Color.h"
#include "Obj2.h"
#include "LinearAlgebra.h"
#include <iostream>
#include <glm/glm.hpp>
#include <vector>
#include <algorithm>

namespace vmr {

	template <typename U>
	struct GradType
	{
		typedef int GType;
	};

	template<>
	struct GradType<float>
	{
		typedef Vector GType;
	};

	template<>
	struct GradType<Vector>
	{
		typedef Matrix GType;
	};


	template < typename U > class Field
	{
	public:

		Field() {}

		virtual ~Field() {}

		typedef U fieldDataType;
		typedef typename GradType<U>::GType fieldGradType;

		virtual const fieldDataType eval(const Vector& P) const { fieldDataType base{}; return base; }
		virtual const fieldGradType grad(const Vector& P) const { fieldGradType base{}; return base; }

	};

	typedef Field<float>* FieldFloatPtr;
	typedef Field<Color>* FieldColorPtr;
	typedef Field<Vector>* FieldVectorPtr;
	typedef Field<Matrix>* FieldMatrixPtr;


	

	//***************************************************************************************************************
	// Scalar Algebra
	//***************************************************************************************************************
	class ScalarFieldAdd : public Field<float>
	{
	public:
		ScalarFieldAdd(FieldFloatPtr f, FieldFloatPtr g)
		{
			e1 = f;
			e2 = g;
		}

		const float eval(Vector& P) const
		{
			std::cout << "\n evallll";
			std::cout << "\n evallll as;ldkfj ;alsk df";
			return e1->eval(P) + e2->eval(P);
		}

	private:
		FieldFloatPtr e1;
		FieldFloatPtr e2;
	};

	inline const ScalarFieldAdd operator+ (Field<float> f, Field<float> g) 
	{
		return ScalarFieldAdd(&f, &g);
	}

	class ScalarFieldAddFloat : public Field<float>
	{
	public:
		ScalarFieldAddFloat(Field<float> * f, float g)
		{
			e1 = f;
			e2 = g;
		}

		float eval(Vector& P)
		{
			return e1->eval(P)+e2;
		}

	private:
		FieldFloatPtr e1;
		float e2;
	};

	class ScalarFieldSubtract : public Field<float>
	{
	public:
		ScalarFieldSubtract(Field<float> f, Field<float> g)
		{
			e1 = f;
			e2 = g;
		}

		float eval(Vector& P)
		{
			return e1.eval(P) - e2.eval(P);
		}

	private:
		Field<float> e1;
		Field<float> e2;
	};

	inline const ScalarFieldSubtract operator- (Field<float> f, Field<float> g)
	{
		return ScalarFieldSubtract(f, g);
	}

	class ScalarFieldMultiply : public Field<float>
	{
	public:
		ScalarFieldMultiply(Field<float> f, Field<float> g)
		{
			e1 = f;
			e2 = g;
		}

		float eval(Vector& P)
		{
			return e1.eval(P) * e2.eval(P);
		}

	private:
		Field<float> e1;
		Field<float> e2;
	};

	inline const ScalarFieldMultiply operator* (Field<float> f, Field<float> g)
	{
		return ScalarFieldMultiply(f, g);
	}

	class ScalarFieldDivide : public Field<float>
	{
	public:
		ScalarFieldDivide(Field<float> f, Field<float> g)
		{
			e1 = f;
			e2 = g;
		}

		float eval(Vector& P)
		{
			return e1.eval(P) - e2.eval(P);
		}

	private:
		Field<float> e1;
		Field<float> e2;
	};

	inline const ScalarFieldDivide operator/ (Field<float> f, Field<float> g)
	{
		return ScalarFieldDivide(f, g);
	}

	class ScalarFieldUnaryMinus : public Field<float>
	{
	public:
		ScalarFieldUnaryMinus(Field<float> f)
		{
			e1 = f;
		}

		float eval(Vector& P)
		{
			return -e1.eval(P);
		}

	private:
		Field<float> e1;
	};

	inline const ScalarFieldUnaryMinus operator- (Field<float> f)
	{
		return ScalarFieldUnaryMinus(f);
	}

	class ScalarFieldInverse : public Field<float>
	{
	public:
		ScalarFieldInverse(Field<float> f)
		{
			e1 = f;
		}

		float eval(Vector& P)
		{
			return 1/e1.eval(P);
		}

	private:
		Field<float> e1;
	};

	inline const ScalarFieldInverse operator/ (int t, Field<float> f)
	{
		return ScalarFieldInverse(f);
	}

	//***************************************************************************************************************
	// Vector Algebra
	//***************************************************************************************************************
	class VectorFieldAdd : public Field<Vector>
	{
	public:
		VectorFieldAdd(Field<Vector> f, Field<Vector> g)
		{
			e1 = f;
			e2 = g;
		}

		Vector eval(Vector& P)
		{
			return e1.eval(P) + e2.eval(P);
		}

	private:
		Field<Vector> e1;
		Field<Vector> e2;
	};

	inline const VectorFieldAdd operator+ (Field<Vector> f, Field<Vector> g)
	{
		return VectorFieldAdd(f, g);
	}

	class VectorFieldSubtract : public Field<Vector>
	{
	public:
		VectorFieldSubtract(Field<Vector> f, Field<Vector> g)
		{
			e1 = f;
			e2 = g;
		}

		Vector eval(Vector& P)
		{
			return e1.eval(P) - e2.eval(P);
		}

	private:
		Field<Vector> e1;
		Field<Vector> e2;
	};

	inline const VectorFieldSubtract operator- (Field<Vector> f, Field<Vector> g)
	{
		return VectorFieldSubtract(f, g);
	}

	class VectorFieldCrossProduct : public Field<Vector>
	{
	public:
		VectorFieldCrossProduct(Field<Vector> f, Field<Vector> g)
		{
			e1 = f;
			e2 = g;
		}

		Vector eval(Vector& P)
		{
			return e1.eval(P) ^ e2.eval(P);
		}

	private:
		Field<Vector> e1;
		Field<Vector> e2;
	};

	inline const VectorFieldCrossProduct operatorx (Field<Vector> f, Field<Vector> g)
	{
		return VectorFieldCrossProduct(f, g);
	}

	class VectorFieldDotProduct : public Field<Vector>
	{
	public:
		VectorFieldDotProduct(Field<Vector> f, Field<Vector> g)
		{
			e1 = f;
			e2 = g;
		}

		float eval(Vector& P)
		{
			return e1.eval(P) * e2.eval(P);
		}

	private:
		Field<Vector> e1;
		Field<Vector> e2;
	};

	inline const VectorFieldDotProduct operator* (Field<Vector> f, Field<Vector> g)
	{
		return VectorFieldDotProduct(f, g);
	}

	class VectorFieldOuterProduct : public Field<Vector> 
	{
	public:
		VectorFieldOuterProduct(Field<Vector> f, Field<Vector> g)
		{
			e1 = f;
			e2 = g;
		}

		Matrix eval(Vector& P)
		{
			return e1.eval(P) & e2.eval(P);
		}

	private:
		Field<Vector> e1;
		Field<Vector> e2;
	};

	inline const VectorFieldOuterProduct operator& (Field<Vector> f, Field<Vector> g)
	{
		return VectorFieldOuterProduct(f, g);
	}

	class VectorFieldScalarMultiplication : public Field<Vector>
	{
	public:
		VectorFieldScalarMultiplication(Field<Vector> f, Field<float> g)
		{
			e1 = f;
			e2 = g;
		}

		Vector eval(Vector& P)
		{
			return e1.eval(P) * e2.eval(P);
		}

	private:
		Field<Vector> e1;
		Field<float> e2;
	};

	inline const VectorFieldScalarMultiplication operator* (Field<Vector> f, Field<float> g)
	{
		return VectorFieldScalarMultiplication(f, g);
	}

	class VectorFieldScalarDivision : public Field<Vector>
	{
	public:
		VectorFieldScalarDivision(Field<Vector> f, Field<float> g)
		{
			e1 = f;
			e2 = g;
		}

		Vector eval(Vector& P)
		{
			return e1.eval(P) / e2.eval(P);
		}

	private:
		Field<Vector> e1;
		Field<float> e2;
	};

	inline const VectorFieldScalarDivision operator/ (Field<Vector> f, Field<float> g)
	{
		return VectorFieldScalarDivision(f, g);
	}

	//***************************************************************************************************************
	// Matrix Algebra
	//***************************************************************************************************************
	class MatrixFieldAddition : public Field<Matrix>
	{
	public:
		MatrixFieldAddition(Field<Matrix> f, Field<Matrix> g)
		{
			e1 = f;
			e2 = g;
		}

		Matrix eval(Vector& P)
		{
			return e1.eval(P) + e2.eval(P);
		}

	private:
		Field<Matrix> e1;
		Field<Matrix> e2;
	};

	inline const MatrixFieldAddition operator+ (Field<Matrix> f, Field<Matrix> g)
	{
		return MatrixFieldAddition(f, g);
	}

	class MatrixFieldSubtraction : public Field<Matrix>
	{
	public:
		MatrixFieldSubtraction(Field<Matrix> f, Field<Matrix> g)
		{
			e1 = f;
			e2 = g;
		}

		Matrix eval(Vector& P)
		{
			return e1.eval(P) - e2.eval(P);
		}

	private:
		Field<Matrix> e1;
		Field<Matrix> e2;
	};

	inline const MatrixFieldSubtraction operator- (Field<Matrix> f, Field<Matrix> g)
	{
		return MatrixFieldSubtraction(f, g);
	}

	class MatrixFieldUnaryMinus : public Field<Matrix>
	{
	public:
		MatrixFieldUnaryMinus(Field<Matrix> f)
		{
			e1 = f;
		}

		Matrix eval(Vector& P)
		{
			return -e1.eval(P);
		}

	private:
		Field<Matrix> e1;
	};

	inline const MatrixFieldUnaryMinus operator- (Field<Matrix> f)
	{
		return MatrixFieldUnaryMinus(f);
	}

	class MatrixFieldMultiplication : public Field<Matrix>
	{
	public:
		MatrixFieldMultiplication(Field<Matrix> f, Field<Matrix> g)
		{
			e1 = f;
			e2 = g;
		}

		Matrix eval(Vector& P)
		{
			return e1.eval(P) * e2.eval(P);
		}

	private:
		Field<Matrix> e1;
		Field<Matrix> e2;
	};

	inline const MatrixFieldMultiplication operator* (Field<Matrix> f, Field<Matrix> g)
	{
		return MatrixFieldMultiplication(f, g);
	}

	class MatrixFieldMVMultiplication : public Field<Matrix>
	{
	public:
		MatrixFieldMVMultiplication(Field<Matrix> f, Field<Vector> g)
		{
			e1 = f;
			e2 = g;
		}

		Vector eval(Vector& P)
		{
			return e1.eval(P) * e2.eval(P);
		}

	private:
		Field<Matrix> e1;
		Field<Vector> e2;
	};

	inline const MatrixFieldMVMultiplication operator* (Field<Matrix> f, Field<Vector> g)
	{
		return MatrixFieldMVMultiplication(f, g);
	}

	class MatrixFieldVMMultiplication : public Field<Matrix>
	{
	public:
		MatrixFieldVMMultiplication(Field<Vector> f, Field<Matrix> g)
		{
			e1 = f;
			e2 = g;
		}

		Vector eval(Vector& P)
		{
			return e1.eval(P) * e2.eval(P);
		}

	private:
		Field<Vector> e1;
		Field<Matrix> e2;
	};

	inline const MatrixFieldVMMultiplication operator* (Field<Vector> f, Field<Matrix> g)
	{
		return MatrixFieldVMMultiplication(f, g);
	}

	/*class MatrixFieldVectorSummation1 : public Field<Matrix> //////////////////////////////////////////////// FINISH
	{
	public:
		MatrixFieldMVMultiplication(Field<Matrix> f, Field<Vector> g)
		{
			e1 = f;
			e2 = g;
		}

		Vector eval(Vector& P)
		{
			// return
		}

	private:
		Field<Matrix> e1;
		Field<Vector> e2;
	};
	
	class MatrixFieldSummation2 : public Field<Matrix>
	{
	public:
		MatrixFieldMVMultiplication(Field<Matrix> f, Field<Vector> g)
		{
			e1 = f;
			e2 = g;
		}

		Vector eval(Vector& P)
		{
			return e1.eval(P) * e2.eval(P);
		}

	private:
		Field<Matrix> e1;
		Field<Vector> e2;
	};
	
	*/

	//***************************************************************************************************************
	// Color Algebra
	//***************************************************************************************************************
	class ColorFieldAddition : public Field<Color>
	{
	public:
		ColorFieldAddition(Field<Color> * f, Field<Color> * g)
		{
			e1 = f;
			e2 = g;
		}

		const Color eval(Vector& P) const
		{
			return e1->eval(P) + e2->eval(P);
		}

	private:
		FieldColorPtr e1;
		FieldColorPtr e2;
	};

	inline const ColorFieldAddition operator+ (Field<Color> &f, Field<Color> &g)
	{
		return ColorFieldAddition(&f, &g);
	}

	class ColorFieldSubtraction : public Field<Color>
	{
	public:
		ColorFieldSubtraction(Field<Color> f, Field<Color> g)
		{
			e1 = f;
			e2 = g;
		}

		Color eval(Vector& P)
		{
			return e1.eval(P) - e2.eval(P);
		}

	private:
		Field<Color> e1;
		Field<Color> e2;
	};

	inline const ColorFieldSubtraction operator- (Field<Color> f, Field<Color> g)
	{
		return ColorFieldSubtraction(f, g);
	}

	class ColorFieldMultiplication : public Field<Color>
	{
	public:
		ColorFieldMultiplication(Field<Color> f, Field<Color> g)
		{
			e1 = f;
			e2 = g;
		}

		Color eval(Vector& P)
		{
			return e1.eval(P) * e2.eval(P);
		}

	private:
		Field<Color> e1;
		Field<Color> e2;
	};

	inline const ColorFieldMultiplication operator* (Field<Color> f, Field<Color> g)
	{
		return ColorFieldMultiplication(f, g);
	}

	class ColorFieldScalarDivision : public Field<Color>
	{
	public:
		ColorFieldScalarDivision(Field<Color> f, Field<float> g)
		{
			e1 = f;
			e2 = g;
		}

		Color eval(Vector& P)
		{
			return e1.eval(P) / e2.eval(P);
		}

	private:
		Field<Color> e1;
		Field<float> e2;
	};

	inline const ColorFieldScalarDivision operator/ (Field<Color> f, Field<float> g)
	{
		return ColorFieldScalarDivision(f, g);
	}

	class ColorFieldScalarMultiplication : public Field<Color>
	{
	public:
		ColorFieldScalarMultiplication(Field<Color> f, Field<float> g)
		{
			e1 = f;
			e2 = g;
		}

		Color eval(Vector& P)
		{
			return e1.eval(P) * e2.eval(P);
		}

	private:
		Field<Color> e1;
		Field<float> e2;
	};

	inline const ColorFieldScalarMultiplication operator* (Field<Color> f, Field<float> g)
	{
		return ColorFieldScalarMultiplication(f, g);
	}

	//***************************************************************************************************************
	// Color Field
	//***************************************************************************************************************
	class ColorField : public Field<Color>
	{

	public:

		ColorField(Color color)
		{
			c = color;
		}

		const Color eval(const Vector& P) const
		{
			return c;
		}

	private:

		Color c;

	};

	//***************************************************************************************************************
	// Gridded Fields
	//***************************************************************************************************************

	float clamp(float v, float low, float high)
	{
		if (v < low)
			return low;
		if (v > high)
			return high;
	}

	bool isIn(Vector pos, Vector boxDim, Vector boxCent)
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

	float clip(float n, float lower, float upper) {
		return std::max(lower, std::min(n, upper));
	}

	class griddedFloatField : public Field<float>
	{

	public:

		griddedFloatField(float * grid, Vector llc, Vector dim,float step)
		{
			g = grid;
			l = llc;
			s = step;
			n = dim;
		}

		const float eval(const Vector& P) const
		{
			//interpolate position at grid
			Vector closestPoint = Vector(floor(((P.X() - l.X()) / s)), floor(((P.Y() - l.Y()) / s)), floor(((P.Z() - l.Z()) / s)));

			float i = (P[0] - l[0]) / s;
			float j = (P[1] - l[1]) / s;
			float k = (P[2] - l[2]) / s;

			float wi = i - (int)(i);
			float wj = j - (int)(j);
			float wk = k - (int)(k);

			i = (int)(i);
			j = (int)(j);
			k = (int)(k);

			int width = ceil(n.X() / s);
			//int width = 30, height = 30, depth = 30;
			int height = ceil(n.Y() / s);
			int depth = ceil(n.Z() / s);

			if (isIn(P, n, (l + n / 2)))
			{
				if (i < width-1 && j < height-1 && k < depth-1)
				{

					int c000 = (i + width * (j + depth * k));
					int c100 = (i + 1) + width * (j + depth * k);
					int c001 = (i)+width * ((j + 1) + depth * k);
					int c010 = (i)+width * (j + depth * (k + 1));
					int c101 = (i + 1) + width * ((j + 1) + depth * (k));
					int c110 = (i + 1) + width * (j + depth * (k + 1));
					int c011 = (i)+width * ((j + 1) + depth * (k + 1));
					int c111 = (i + 1) + width * ((j + 1) + depth * (k + 1));

					float gx = g[c000] * (1.0 - wi) * (1.0 - wj) * (1.0 - wk)
						+ g[c100] * wi * (1.0 - wj) * (1.0 - wk)
						+ g[c001] * (1.0 - wi) * wj * (1.0 - wk)
						+ g[c010] * (1.0 - wi) * (1.0 - wj) * (wk)
						+g[c101] * (wi) * (wj) * (1.0 - wk)
						+ g[c110] * (wi) * (1.0 - wj) * (wk)
						+g[c011] * (1.0 - wi) * (wj) * (wk)
						+g[c111] * wi * wj * wk;
					return gx;
				}
				else
				{
					return 0;
				}
			}
			else
			{
				return 0;
			}
		
		}

		const Vector grad(const Vector&P)const 
		{
			float tmpx = eval(P + Vector(1, 0, 0)*0.001) - eval(P-Vector(1,0,0)*0.001) / (2 * 0.001);
			float tmpy = eval(P + Vector(0, 1, 0)*0.001) - eval(P-Vector(0,1,0)*0.001) / (2 * 0.001);
			float tmpz = eval(P + Vector(0, 0, 1)*0.001) - eval(P-Vector(0,0,1)*0.001) / (2 * 0.001);
			return Vector(tmpx, tmpy, tmpz);
		}

	private:

		float * g;
		Vector l;
		float s;
		Vector n;

	};

	class griddedVectorField : public Field<Vector>
	{

	public:

		griddedVectorField(float * grid, Vector llc, Vector dim, float step)
		{
			g = grid;
			l = llc;
			s = step;
			n = dim;
		}

		const Vector eval(const Vector& P) const
		{
			//interpolate position at grid
			Vector closestPoint = Vector(floor(((P.X() - l.X()) / s)), floor(((P.Y() - l.Y()) / s)), floor(((P.Z() - l.Z()) / s)));

			float i = (int)(P[0] - l[0]) / s;
			float j = (int)(P[1] - l[1]) / s;
			float k = (int)(P[2] - l[2]) / s;

			int w = ceil(n.X() / s);
			int height = ceil(n.Y() / s);
			int depth = ceil(n.Z() / s);

			if (isIn(P, n, (l + n / 2)))
			{
				if (i < w - 1 && j < height - 1 && k < depth - 1)
				{

					int c0 = i + (j * w) + (k * w * w) + (0 * w * w * w);
					int c1 = i + (j * w) + (k * w * w) + (1 * w * w * w);
					int c2 = i + (j * w) + (k * w * w) + (2 * w * w * w);
					return Vector(g[c0], g[c1], g[c2]);
				}
				else
				{
					return Vector(0,0,0);
				}
			}
			else
			{
				return Vector(0,0,0);
			}

		}

	private:

		float * g;
		Vector l;
		float s;
		Vector n;

	};

	class LevelSetField : public Field<float> {

	public:
		LevelSetField(std::vector<Face> faces) {
			myFaces = faces;
		}

		const float eval(const Vector&P)const {

			Vector p0, p1, p2; // vertices
			Vector e1, e2, e3; // edges
			Vector vn0, vn1, vn2; //vertex normals
			Vector pN; //plane normal
			Vector r; 

			float minDist = -1000;
			for (int i = 0; i < myFaces.size(); i++) {

				p0 = myFaces.at(i).getVertex0();
				p1 = myFaces.at(i).getVertex1();
				p2 = myFaces.at(i).getVertex2();

				//vn0 = Vector(norms[0][i][0], norms[0][i][1], norms[0][i][2]);
				//vn1 = Vector(norms[0][i + 1][0], norms[0][i + 1][1], norms[0][i + 1][2]);
				//vn2 = Vector(norms[0][i + 2][0], norms[0][i + 2][1], norms[0][i + 2][2]);

				e1 = (p1 - p0);
				e2 = (p2 - p0);
				e3 = (p2 - p1);

				pN = (e1^e2).unitvector();
				float u, v;
				r = P - p0 - (pN * (pN*(P - p0)));
				v = ((e1^r)*(e1^e2)) / ((e1^e2).magnitude()*(e1^e2).magnitude());
				u = ((e2 ^ r)*(e2^e1)) / ((e2^e1).magnitude()*(e2^e1).magnitude());
				
				//faces
				float dist;
				float det;
				if (u >= 0 && u <= 1 && v >= 0 && v <= 1 && (u + v) >= 0 && (u + v) <= 1) 
				{
					dist = pN * (p0 - P);
					if (abs(dist) < abs(minDist))
						minDist = dist;
				}
				
				//edge 1
				float q = (e1*(P - p0)) / (e1.magnitude()*e1.magnitude());
				if (q >= 0 && q <= 1) 
				{
					dist = ((p0 + (q*e1)) - P).magnitude();
					det = (P - (p0 + (q*e1)))*pN;
					if (det > 0)
						dist = -dist;
					if (abs(dist) < abs(minDist))
						minDist = dist;
				}

				//edge2
				q = (e2*(P - p0)) / (e2.magnitude()*e2.magnitude());
				if (q >= 0 && q <= 1) {
					dist = ((p0 + (q*e2)) - P).magnitude();
					det = (P - (p0 + (q*e2)))*pN;
					if (det > 0)
						dist = -dist;
					if (abs(dist) < abs(minDist))
						minDist = dist;
				}

				//edge3
				q = (e3*(P - p1)) / (e3.magnitude()*e3.magnitude());
				if (q >= 0 && q <= 1) {
					dist = ((p1 + (q*e3)) - P).magnitude();
					det = (P - (p1 + (q*e3)))*pN;
					if (det > 0)
						dist = -dist;
					if (abs(dist) < abs(minDist))
						minDist = dist;
				}

				//p0
				dist = (P - p0).magnitude();
				det = (P - p0)*pN;
				if (det > 0) 
					dist = -dist;
				if (abs(dist) < abs(minDist))
					minDist = dist;
				
				//p1
				dist = (P - p1).magnitude();
				det = (P - p1)*pN;
				if (det > 0) 
					dist = -dist;
				if (abs(dist) < abs(minDist))
					minDist = dist;

				//p2
				dist = (P - p2).magnitude();
				det = (P - p2)*pN;
				if (det > 0) 
					dist = -dist;
				if (abs(dist) < abs(minDist))
					minDist = dist;

			}

			return minDist;
		}

	private:
		std::vector<Face> myFaces;
	};

	class LevelSetFieldFaces : public Field<float> {

	public:
		LevelSetFieldFaces(std::vector<Face> faces) {
			myFaces = faces;
		}

		const float eval(const Vector&P)const {

			Vector p0, p1, p2;
			Vector edge1, edge2, edge3;
			Vector vnorm0, vnorm1, vnorm2;
			float u, v;
			Vector planeNorm;
			Vector r;

			float minDist = -100000;
			bool hasHit = false;
			//std::cout << "in eval";
			for (int i = 0; i < myFaces.size(); i++) {

				float d;
				p0 = myFaces.at(i).getVertex0();
				p1 = myFaces.at(i).getVertex1();
				p2 = myFaces.at(i).getVertex2();

				edge1 = (p1 - p0);
				edge2 = (p2 - p0);
				edge3 = (p2 - p1);

				vnorm0 = myFaces.at(i).getNorm0();
				vnorm1 = myFaces.at(i).getNorm1();
				vnorm2 = myFaces.at(i).getNorm2();

				planeNorm = (edge1^edge2).unitvector();

				r = P - p0 - (planeNorm * (planeNorm*(P - p0)));
				v = ((edge1^r)*(edge1^edge2)) / ((edge1^edge2).magnitude()*(edge1^edge2).magnitude());
				u = ((edge2 ^ r)*(edge2^edge1)) / ((edge2^edge1).magnitude()*(edge2^edge1).magnitude());

				//check face
				if (u >= 0 && u <= 1 && v >= 0 && v <= 1 && (u + v) >= 0 && (u + v) <= 1) {
					d = planeNorm * (p0 - P);
					if (abs(d) < abs(minDist))
						minDist = d;
				}

				//std::cout << "\n mindist: " << minDist;
				//check edges
				float q = (edge1*(P - p0)) / (edge1.magnitude()*edge1.magnitude());
				if (q >= 0 && q <= 1) {
					d = ((p0 + (q*edge1)) - P).magnitude();
					if ((P - (p0 + (q*edge1)))*planeNorm > 0)
						d = -d;
					if (abs(d) < abs(minDist))
						minDist = d;
				}

				q = (edge2*(P - p0)) / (edge2.magnitude()*edge2.magnitude());
				if (q >= 0 && q <= 1) {
					d = ((p0 + (q*edge2)) - P).magnitude();
					if ((P - (p0 + (q*edge2)))*planeNorm > 0)
						d = -d;
					if (abs(d) < abs(minDist))
						minDist = d;
				}


				q = (edge3*(P - p1)) / (edge3.magnitude()*edge3.magnitude());
				if (q >= 0 && q <= 1) {
					d = ((p1 + (q*edge3)) - P).magnitude();
					if ((P - (p1 + (q*edge3)))*planeNorm > 0)
						d = -d;
					if (abs(d) < abs(minDist))
						minDist = d;
				}

				float h, h1, h2, h3;
				//check vertices
				h1 = (P - p0).magnitude();
				if ((P - p0)*planeNorm > 0) //get sign
					h1 = -h1;

				h2 = (P - p1).magnitude();
				if ((P - p1)*planeNorm > 0) //get sign
					h2 = -h2;

				h3 = (P - p2).magnitude();
				if ((P - p2)*planeNorm) //get sign
					h3 = -h3;

				h = fmin(abs(h1), (fmin(abs(h2), abs(h3))));

				if (abs(h) < minDist)
					h = minDist;

			}

			return minDist;
		}

		//const int grad(const Vector&P)const {
		//    return 0;
		//}
		const float& getsize() const { return myFaces.size(); }

	private:
		std::vector<Face> myFaces;
	};

	class LevelSetGradField : public Field<Vector> {

	public:
		LevelSetGradField(std::vector < glm::vec3 > * vertices, std::vector<glm::vec3> * normals) {
			verts = vertices;
			norms = normals;
		}

		const Vector eval(const Vector&P)const {

			Vector p0, p1, p2; // vertices
			Vector e1, e2, e3; // edges
			Vector vn0, vn1, vn2; //vertex normals
			Vector pN; //plane normal
			Vector r;
			Vector iN = Vector(0, 0, 0);

			float minDist = -1000;
			for (int i = 0; i < verts->size(); i += 3)
			{

				p0 = Vector(verts[0][i][0], verts[0][i][1], verts[0][i][2]);
				p1 = Vector(verts[0][i + 1][0], verts[0][i + 1][1], verts[0][i + 1][2]);
				p2 = Vector(verts[0][i + 2][0], verts[0][i + 2][1], verts[0][i + 2][2]);

				vn0 = Vector(norms[0][i][0], norms[0][i][1], norms[0][i][2]);
				vn1 = Vector(norms[0][i + 1][0], norms[0][i + 1][1], norms[0][i + 1][2]);
				vn2 = Vector(norms[0][i + 2][0], norms[0][i + 2][1], norms[0][i + 2][2]);

				e1 = (p1 - p0);
				e2 = (p2 - p0);
				e3 = (p2 - p1);

				pN = (e1^e2).unitvector();
				float u, v;
				r = P - p0 - (pN * (pN*(P - p0)));
				v = ((e1^r)*(e1^e2)) / ((e1^e2).magnitude()*(e1^e2).magnitude());
				u = ((e2 ^ r)*(e2^e1)) / ((e2^e1).magnitude()*(e2^e1).magnitude());

				//faces
				float dist;
				float det;
				if (u >= 0 && u <= 1 && v >= 0 && v <= 1 && (u + v) >= 0 && (u + v) <= 1)
				{
					dist = pN * (p0 - P);
					if (abs(dist) < abs(minDist))
					{
						minDist = dist;
						iN = -pN;
					}
				}

				//edge 1
				float q = (e1*(P - p0)) / (e1.magnitude()*e1.magnitude());
				if (q >= 0 && q <= 1)
				{
					dist = ((p0 + (q*e1)) - P).magnitude();
					det = (P - (p0 + (q*e1)))*pN;
					if (det > 0)
						dist = -dist;
					if (abs(dist) < abs(minDist))
					{
						minDist = dist;
						iN = -pN;
					}
				}

				//edge2
				q = (e2*(P - p0)) / (e2.magnitude()*e2.magnitude());
				if (q >= 0 && q <= 1) {
					dist = ((p0 + (q*e2)) - P).magnitude();
					det = (P - (p0 + (q*e2)))*pN;
					if (det > 0)
						dist = -dist;
					if (abs(dist) < abs(minDist))
					{
						minDist = dist;
						iN = -pN;
					}
				}

				//edge3
				q = (e3*(P - p1)) / (e3.magnitude()*e3.magnitude());
				if (q >= 0 && q <= 1) {
					dist = ((p1 + (q*e3)) - P).magnitude();
					det = (P - (p1 + (q*e3)))*pN;
					if (det > 0)
						dist = -dist;
					if (abs(dist) < abs(minDist))
					{
						minDist = dist;
						iN = -pN;
					}
				}

				//p0
				dist = (P - p0).magnitude();
				det = (P - p0)*pN;
				if (det > 0)
					dist = -dist;
				if (abs(dist) < abs(minDist))
				{
					minDist = dist;
					iN = -pN;
				}

				//p1
				dist = (P - p1).magnitude();
				det = (P - p1)*pN;
				if (det > 0)
					dist = -dist;
				if (abs(dist) < abs(minDist))
				{
					minDist = dist;
					iN = -pN;
				}

				//p2
				dist = (P - p2).magnitude();
				det = (P - p2)*pN;
				if (det > 0)
					dist = -dist;
				if (abs(dist) < abs(minDist))
				{
					minDist = dist;
					iN = -pN;
				}

			}

			return iN;
			//return minDist;
		}

	private:
		std::vector<glm::vec3> *verts;
		std::vector<glm::vec3> *norms;
	};


	class GradVectorField : public Field<Vector>
	{

	public:

		GradVectorField(FieldFloatPtr f1)
		{
			f = f1;
		}

		const Vector eval(const Vector& P) const
		{
			return f->grad(P);
		}

	private:
		FieldFloatPtr f;
	};

	class FSPNVectorField : public Field<Vector>
	{

	public:

		FSPNVectorField(FieldFloatPtr f1,FieldFloatPtr noise,float delta)
		{
			f = f1;
			n = noise;
			d = delta;
			de = Vector(d, d, d);
		}

		const Vector eval(const Vector& P) const
		{
			return Vector(n->eval(P), n->eval(P + de), n->eval(P - de));
		}

	private:
		FieldFloatPtr f;
		FieldFloatPtr n;
		float d;
		Vector de;
	};

	class FloatAdvect : public Field<float>
	{

	public:

		FloatAdvect(FieldFloatPtr f1,FieldVectorPtr v1,float deltaT)
		{
			f = f1;
			v = v1;
			d = deltaT;
		}

		const float eval(const Vector& P) const
		{
			//std::cout << "\n" << v->eval(P)[0] << " " << v->eval(P)[1] << " " << v->eval(P)[2];
			Vector pp = P - v->eval(P)*d;
			return f->eval(pp);
		}

	private:
		FieldFloatPtr f;
		FieldVectorPtr v;
		float d;
	};

	class FloatAdvectItr : public Field<float>
	{

	public:

		FloatAdvectItr(FieldFloatPtr f1, FieldVectorPtr v1, float deltaT,int iter)
		{
			f = f1;
			v = v1;
			d = deltaT;
			it = iter;
		}

		const float eval(const Vector& P) const
		{
			
			Vector pp = P - v->eval(P)*d;

			for (int i = 0; i < it; i++)
			{
				pp = pp - (v->eval(pp)*d);
			}

			return f->eval(pp);
		}

	private:
		FieldFloatPtr f;
		FieldVectorPtr v;
		float d;
		int it;
	};



}

#endif //FIELD_H
