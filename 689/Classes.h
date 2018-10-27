#pragma once

#ifndef VMR_CLASSES_H
#define VMR_CLASSES_H

#include "Vector.h"
#include "Color.h"

namespace vmr
{

	class Light
	{

	public:

		Light(Vector location, Color color)
		{
			l = location;
			c = color;
		}

		Vector location()
		{
			return l;
		}

		Color color()
		{
			return c;
		}

	private:

		Vector l;
		Color c;

	};

}
#endif // !VMR_CLASSES_H
