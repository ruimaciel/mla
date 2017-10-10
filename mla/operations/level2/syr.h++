#ifndef MLA_OPERATIONS_LEVEL2_SYR_HPP
#define MLA_OPERATIONS_LEVEL2_SYR_HPP

#include <type_traits>

#include <mla/LAException.h++>

#include <mla/matrix/all.h++>
#include <mla/vector/all.h++>

#include <mla/MatrixCursor.h++>

namespace mla {

/**
 * Implements the level 2 BLAS functions
 * {A} := A + alpha*{x}{x}^T
 *
 * http://www.netlib.org/blas/#_reference_blas_version_3_5_0
 **/

template<typename Scalar, template<typename> class VectorPolicyX, template<typename> class MatrixPolicyA>
void 
syr(Scalar alpha, VectorPolicyX<Scalar> &x, MatrixPolicyA<Scalar> &A)
{
	if( !A.isSquare()) 
	{
		throw LAException("level2::syr: A must be a square matrix");
	}
	if( x.size() != A.columns() )
	{
		throw LAException("level2::syr: incompatible sizes between A and x");
	}

	Scalar temp_value;
	auto cursor_xi = x.cursor();

	while( !cursor_xi.at_end() )
	{
		auto cursor_xj = cursor_xi;
		size_t i = cursor_xi.current();
		size_t j = cursor_xj.current();

		// first, the diagonal elements
		temp_value = alpha*cursor_xi.element()*cursor_xj.element();
		A(i,j) += temp_value;

		cursor_xj.next();

		// next, the off-diagonal elements
		while( !cursor_xj.at_end() )
		{
			j = cursor_xj.current();

			temp_value = alpha*cursor_xi.element()*cursor_xj.element();
			A(i,j) += temp_value;
			A(j,i) += temp_value;

			cursor_xj.next();
		}

		cursor_xi.next();
	}

}



}	// namespace mla

#endif
