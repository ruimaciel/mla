#ifndef MLA_OPERATIONS_LEVEL3_SYRK_HPP
#define MLA_OPERATIONS_LEVEL3_SYRK_HPP

#include <type_traits>
#include <mla/matrix/all.h++>
#include <mla/vector/all.h++>

namespace mla {

/**
 * Implements the level 3 BLAS functions SYRK
 * {C} := alpha*[A]*[A]^T + beta*[C]
 *
 * http://www.netlib.org/blas/#_level_3
 **/


template<typename Scalar, template<typename> class MatrixAPolicy, template<typename> class MatrixCPolicy>
void
syrk(Scalar const alpha, MatrixAPolicy<Scalar> const &A, Scalar const beta, MatrixCPolicy<Scalar> &C)
{
	static_assert(true, "Generic level3::gemm not implemented yet");
}


template<typename Scalar>
void
syrk(Scalar const alpha, matrix::DenseRowMajor<Scalar> const &A, Scalar const beta, matrix::DenseRowMajor<Scalar> &C)
{
	if(C.rows() != C.columns())
	{
		throw LAException("syrk: C isn't square");
	}
	if(A.rows() != C.rows())
	{
		throw LAException("syrk: A.rows() != C.rows()");
	}

	// proceed 
	for(size_t i = 0; i < A.rows(); i++)
	{
		Scalar value = 0.0f;

		// update diagonal elements
		for(size_t k = 0; k < A.columns(); k++)
		{
			value += A.getValue(i,k)*A.getValue(i,k);
		}
		C(i,i)  = alpha*value + beta*C(i,i);

		// update off-diagonal elements
		for(size_t j = i+1; j < C.columns(); j++)
		{
			// C.setValue(i,j, (Scalar) 0.0f);

			value = 0.0f;
			for(size_t k = 0; k < A.columns(); k++)
			{
				value += A.getValue(i,k)*A.getValue(j,k);
			}

			C(i,j)  = alpha*value + beta*C(i,j);
			C(j,i)  = alpha*value + beta*C(j,i);
		}
	}
}


}	// mla

#endif
