#ifndef MLA_OPERATIONS_LEVEL3_GEMM_HPP
#define MLA_OPERATIONS_LEVEL3_GEMV_HPP

#include <type_traits>
#include <mla/matrix/all.h++>
#include <mla/vector/all.h++>

namespace mla {

/**
 * Implements the level 3 BLAS functions GEMM
 * {C} := a[A][B] + b[C]
 *
 * http://www.netlib.org/blas/#_level_3
 **/


template<typename Scalar, template<typename> class MatrixAPolicy, template<typename> class MatrixBPolicy, template<typename> class MatrixCPolicy>
void
gemm(Scalar const alpha, MatrixAPolicy<Scalar> const &A, MatrixBPolicy<Scalar> const &B, Scalar const beta, MatrixCPolicy<Scalar> &C)
{
	static_assert(true, "Generic level3::gemm not implemented yet");
}


template<typename Scalar>
void
gemm(Scalar const alpha, matrix::DenseRowMajor<Scalar> const &A, matrix::DenseRowMajor<Scalar> const &B, Scalar const beta, matrix::DenseRowMajor<Scalar> &C)
{
	if(A.columns() != B.rows())
	{
		throw LAException("gemm: A.rows() != B.columns()");
	}
	if(C.rows() != A.rows())
	{
		throw LAException("gemm: C.rows() != A.rows()");
	}
	if(C.columns() != B.columns())
	{
		throw LAException("gemm: C.columns() != B.columns()");
	}

	// proceed 
	for(size_t i = 0; i < A.rows(); i++)
	{
		for(size_t j = 0; j < B.columns(); j++)
		{
			// C.setValue(i,j, (Scalar) 0.0f);
			C(i,j) *= beta;

			Scalar value = 0.0f;
			for(size_t k = 0; k < A.columns(); k++)
			{
				value += A.getValue(i,k)*B.getValue(k,j);
			}
			C(i,j)  += alpha*value;
		}
	}
}


}	// mla

#endif
