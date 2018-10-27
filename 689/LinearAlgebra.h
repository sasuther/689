#pragma once 

#ifndef VMR_LINEARALGEBRA_H
#define VMR_LINEARALGEBRA_H

#include "Vector.h"
#include "Matrix.h"

namespace vmr {

	// Matrix times a Vector
	const Vector operator* (const Vector& v, const Matrix& m);
	const Vector operator* (const Matrix& m, const Vector& v);

	// outer product
	const Matrix operator& (const Vector& v1, const Vector& v2);
	// outer product in place
	void outer_product(const Vector& v1, const Vector& v2, Matrix& m);

	const Matrix exp(const Matrix& m);
	const Matrix inverse(const Matrix& m);
	const double det(const Matrix& m);
	const double trace(const Matrix& m);

	// construction a rotation matrix
	const Matrix rotation(const Vector& axis, const double angle);

	const Matrix unitMatrix();

	const static Matrix PauliMatrix[3] = {
									 Matrix(0,0,0,
											 0,0,1,
											 0,-1,0),
									 Matrix(0,0,-1,
											 0,0,0,
											 1,0,0),
									 Matrix(0,1,0,
											 -1,0,0,
											 0,0,0)
	};
}


#endif

