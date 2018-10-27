#pragma once
#ifndef VMR_CSG_H
#define VMR_CSG_H

#include "Field.h"
#include "Vector.h"

namespace vmr
{

	class Union : public Field<float>
	{

	public:

		Union(Field<float> * f1, Field<float> * g1)
		{
			f = f1;
			g = g1;
		}

		const float eval(const Vector& P) const
		{
			return fmax(f->eval(P), g->eval(P));
		}

	private:

		FieldFloatPtr f;
		FieldFloatPtr g;

	};

	class Intersection : public Field<float>
	{

	public:

		Intersection(Field<float> * f1, Field<float> * g1)
		{
			f = f1;
			g = g1;
		}

		const float eval(const Vector& P) const
		{
			return fmin(f->eval(P), g->eval(P));
		}

	private:

		FieldFloatPtr f;
		FieldFloatPtr g;

	};

	class Cutout : public Field<float>
	{

	public:

		Cutout(Field<float> * f1, Field<float> * g1)
		{
			f = f1;
			g = g1;
		}

		const float eval(const Vector& P) const
		{
			return fmin(f->eval(P), -(g->eval(P)));
		}

	private:

		FieldFloatPtr f;
		FieldFloatPtr g;

	};

	class Shell : public Field<float>
	{

	public:
		Shell(FieldFloatPtr f, float g) {
			elem1 = f;
			elem2 = g;
		}

		const float eval(const Vector&P)const {
			return fmin((elem1->eval(P) + elem2), (elem1->eval(P) - elem2));
		}

	private:
		FieldFloatPtr elem1;
		float elem2;

	};


	class ColorIntersection : public Field<Color>
	{

	public:

		ColorIntersection(Field<Color> * f1, Field<float> * g1)
		{
			f = f1;
			g = g1;
		}

		const Color eval(const Vector& P) const
		{
			if (g->eval(P) > 0)
				return f->eval(P);
			else
				return Color(0, 0, 0, 0);
		}

	private:

		FieldColorPtr f;
		FieldFloatPtr g;

	};

	class ColorUnion : public Field<Color>
	{

	public:

		ColorUnion(FieldColorPtr f1, Field<Color> * g1)
		{
			f = f1;
			g = g1;
		}

		const Color eval(Vector& P) const
		{
			//ColorFieldAddition tmp = ColorFieldAddition(f, g);
			//return tmp.eval(P);
			return (f->eval(P) + g->eval(P));
			//return (f+g).eval(P);
		}

	private:

		FieldColorPtr f;
		FieldColorPtr g;

	};

	class ColorTestUnion : public Field<Color>
	{

	public:
		ColorTestUnion(FieldColorPtr f, FieldColorPtr g) {
			elem1 = f;
			elem2 = g;
		}

		const Color eval(const Vector&P)const {
			Color c1 = elem1->eval(P);
			Color c2 = elem2->eval(P);
			return (c1 + c2);
		}

	private:
		FieldColorPtr elem1;
		FieldColorPtr elem2;

	};

}

#endif // !VMR_CSG_H
