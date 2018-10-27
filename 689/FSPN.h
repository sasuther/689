#pragma once
#ifndef VMR_FSPN_H
#define VMR_FSPN_H

#include "PerlinNoise.h"
#include "Vector.h"

namespace vmr
{

	struct FSPN
	{
		float pscale;
		float freq;
		float fJump;
		float rough;
		float oct;
	};

	struct Wisp
	{
		Vector pos;
		float pscale;
		float density; 
		float clump; 
		Vector offset; 
		float offScale;
		int numDots;
	};
}



#endif // !VMR_FSPN_H
