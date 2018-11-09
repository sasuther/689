#pragma once

#ifndef VMR_TRANSFORMS_H
#define VMR_TRANSFORMS_H

#include "Field.h"
#include "Implicit.h"
#include "Vector.h"
#include "Matrix.h"
#include "LinearAlgebra.h"

namespace vmr
{
	//***************************************************************************************************************
	// Translate
	//***************************************************************************************************************
	
	class FloatTranslate : public Field<float>
	{

	public:

		FloatTranslate(Field<float> * field, Vector delta)
		{
			f = field;
			d = delta;
		}

		const float eval(const Vector& P) const
		{
			return f->eval(P-d);
		}

	private:

		FieldFloatPtr f;
		Vector d;

	};

	class VectorTranslate : public Field<Vector>
	{

	public:

		VectorTranslate(Field<Vector> * field, Vector delta)
		{
			f = field;
			d = delta;
		}

		const Vector eval(const Vector& P) const
		{
			return f->eval(P - d);
		}

	private:

		FieldVectorPtr f;
		Vector d;

	};

	class MatrixTranslate : public Field<Matrix>
	{

	public:

		MatrixTranslate(Field<Matrix> * field, Vector delta)
		{
			f = field;
			d = delta;
		}

		const Matrix eval(const Vector& P) const
		{
			return f->eval(P - d);
		}

	private:

		FieldMatrixPtr f;
		Vector d;

	};

	//***************************************************************************************************************
	// Scale
	//***************************************************************************************************************

	class FloatScale : public Field<float>
	{

	public:

		FloatScale(Field<float> * field, Vector delta)
		{
			f = field;
			d = delta;
		}

		const float eval(const Vector& P) const
		{
			return f->eval(Vector(P[0] / d[0], P[1] / d[1], P[2] / d[2]));
		}

	private:

		FieldFloatPtr f;
		Vector d;

	};

	class VectorScale : public Field<Vector>
	{

	public:

		VectorScale(Field<Vector> * field, float delta)
		{
			f = field;
			d = delta;
		}

		const Vector eval(const Vector& P) const
		{
			return f->eval(P / d)*d;
		}

	private:

		FieldVectorPtr f;
		float d;

	};

	class MatrixScale : public Field<Matrix>
	{

	public:

		MatrixScale(Field<Matrix> * field, float delta)
		{
			f = field;
			d = delta;
		}

		const Matrix eval(const Vector& P) const
		{
			return f->eval(P / d)*d*d;
		}

	private:

		FieldMatrixPtr f;
		float d;

	};

	//***************************************************************************************************************
	// Rotation
	//***************************************************************************************************************

	class FloatRotation : public Field<float>
	{

	public:

		FloatRotation(Field<float> * field, Vector axis, float angle)
		{
			f = field;
			a = axis.unitvector();
			an = angle;
		}

		const float eval(const Vector& P) const
		{
			Matrix r = rotation(a, an);
			return f->eval(r.inverse()*P);
		}

	private:

		FieldFloatPtr f;
		Vector a;
		float an;

	};

}

#endif // !VMR_TRANSFORMS_H
