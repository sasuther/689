
#pragma once 

#ifndef FIELD_H
#define FIELD_H

#include "Vector.h"
#include "Matrix.h"
#include "Color.h"
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
		ScalarFieldAdd(Field<float> f, Field<float> g)
		{
			e1 = f;
			e2 = g;
		}

		float eval(Vector& P)
		{
			return e1.eval(P) + e2.eval(P);
		}

	private:
		Field<float> e1;
		Field<float> e2;
	};

	inline const ScalarFieldAdd operator+ (Field<float> f, Field<float> g)
	{
		return ScalarFieldAdd(f, g);
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

			float width = (int)(n.X() / s);
			float height = (int)(n.Y() / s);
			float depth = (int)(n.Z() / s);

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
					return -1;
				}
			}
			else
			{
				return -1;
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
		LevelSetField(std::vector < glm::vec3 > * vertices, std::vector<glm::vec3> * normals) {
			verts = vertices;
			norms = normals;
		}

		const float eval(const Vector&P)const {

			Vector p0, p1, p2; // vertices
			Vector e1, e2, e3; // edges
			Vector vn0, vn1, vn2; //vertex normals
			Vector pN; //plane normal
			Vector r; 

			float minDist = -1000;
			for (int i = 0; i < verts->size(); i+=3) 
			{

				p0 = Vector(verts[0][i][0], verts[0][i][1], verts[0][i][2]);
				p1 = Vector(verts[0][i+1][0], verts[0][i+1][1], verts[0][i+1][2]);
				p2 = Vector(verts[0][i+2][0], verts[0][i+2][1], verts[0][i+2][2]);

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
						minDist = dist;
				}
				
				//edge 1
				float q = (e1*(P - p0)) / (e1.magnitude()*e1.magnitude());
				if (q >= 0 && q <= 1) 
				{
					dist = ((p0 + (q*e1)) - P).magnitude();
					det = (P - (p0 + (q*e1)))*vn0;
					if (det > 0)
						dist = -dist;
					if (abs(dist) < abs(minDist))
						minDist = dist;
				}

				//edge2
				q = (e2*(P - p0)) / (e2.magnitude()*e2.magnitude());
				if (q >= 0 && q <= 1) {
					dist = ((p0 + (q*e2)) - P).magnitude();
					det = (P - (p0 + (q*e2)))*vn0;
					if (det > 0)
						dist = -dist;
					if (abs(dist) < abs(minDist))
						minDist = dist;
				}

				//edge3
				q = (e3*(P - p1)) / (e3.magnitude()*e3.magnitude());
				if (q >= 0 && q <= 1) {
					dist = ((p1 + (q*e3)) - P).magnitude();
					det = (P - (p1 + (q*e3)))*vn1;
					if (det > 0)
						dist = -dist;
					if (abs(dist) < abs(minDist))
						minDist = dist;
				}

				//p0
				dist = (P - p0).magnitude();
				det = (P - p0)*vn0;
				if (det > 0) 
					dist = -dist;
				if (abs(dist) < abs(minDist))
					minDist = dist;
				
				//p1
				dist = (P - p1).magnitude();
				det = (P - p1)*vn1;
				if (det > 0) 
					dist = -dist;
				if (abs(dist) < abs(minDist))
					minDist = dist;

				//p2
				dist = (P - p2).magnitude();
				det = (P - p2)*vn2;
				if (det > 0) 
					dist = -dist;
				if (abs(dist) < abs(minDist))
					minDist = dist;

			}

			return minDist;
		}

	private:
		std::vector<glm::vec3> *verts;
		std::vector<glm::vec3> *norms;
	};


}

#endif //FIELD_H
