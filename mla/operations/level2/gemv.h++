#ifndef MLA_OPERATIONS_LEVEL2_GEMV_HPP
#define MLA_OPERATIONS_LEVEL2_GEMV_HPP

#include <type_traits>

#include <mla/LAException.h++>

#include <mla/matrix/all.h++>
#include <mla/vector/all.h++>

#include <mla/MatrixCursor.h++>

namespace mla {

/**
 * Implements the level 1 BLAS functions
 * {y} := a[A]{x} + b{y}
 *
 * http://www.netlib.org/blas/#_reference_blas_version_3_5_0
 **/

template<typename Scalar, template<typename> class MatrixPolicy, template<typename> class VectorPolicyX, template<typename> class VectorPolicyY>
void 
gemv(Scalar const a, MatrixPolicy<Scalar> &A, VectorPolicyX<Scalar> &x, Scalar const b, VectorPolicyY<Scalar> &y)
{
	if( A.columns() != x.size() )
	{
		throw LAException("level2::gemv: incompatible sizes between A and x");
	}
	if( A.rows() != y.size() )
	{
		throw LAException("level2::gemv: incompatible sizes between A and y");
	}

	typename MatrixPolicy<Scalar>::Cursor A_cursor = A.cursor();
	typename VectorPolicyX<Scalar>::Cursor x_cursor = x.cursor();

	A_cursor.reset();

	while( !A_cursor.at_end_of_rows() )
	{
		x_cursor.reset();

		Scalar Yval = y.getValue(A_cursor.current_row());

		while( !A_cursor.at_end_of_rows() )
		{
			if(A_cursor.current_column()  == x_cursor.current())
			{
				Yval = a*A_cursor.element()*x_cursor.element() + b*Yval;

				x_cursor.next();
				A_cursor.increment_column();

				if( x_cursor.at_end() ||  A_cursor.at_end_of_current_row() )
				{
					y.setValue( A_cursor.current_row(), Yval );
					A_cursor.start_next_row_nn();
					break;
				}
			}
			else if(A_cursor.current_column() > x_cursor.current())
			{
				x_cursor.next();
				if( x_cursor.at_end() )
				{
					y.setValue( A_cursor.current_row(), Yval );
					A_cursor.start_next_column_nn();
					break;
				}
			}
			else
			{
				A_cursor.increment_column();
				if(A_cursor.at_end_of_current_row() )
				{
					y.setValue( A_cursor.current_row(), Yval );
					A_cursor.start_next_row_nn();
					break;
				}
			}
		}

	}
}


/**
The default matrix-vector operator
@param	A	a matrix, instance of class mla::matrix::DenseRowMajor<Scalar>
@param	x	a vector, instance of class VectorStoragePolicy<Scalar>
@param	y	a vector, instance of class VectorStoragePolicy<Scalar>
**/
template<typename Scalar, template<typename> class VectorPolicyX, template<typename> class VectorPolicyY > 
void
gemv(Scalar const a, matrix::DenseRowMajor<Scalar> &A, VectorPolicyX<Scalar> &x, Scalar const b, VectorPolicyY<Scalar> &y)
{
	if(A.columns() != x.size())
	{
		throw LAException("mv: m.columns() != v.size()");
	}
	if(A.rows() != y.size())
	{
		y.resize( A.rows() );
	}

	auto x_cursor = x.cursor();


	for(size_t i = 0; i < A.rows(); i++)
	{
		x_cursor.reset();

		Scalar Yval = y.getValue(i);

		while( !x_cursor.at_end() )
		{
			Yval = a*A.getValue(i, x_cursor.current())*x_cursor.element() + b*Yval;
			x_cursor.next();
		}

		y.setValue( i, Yval);
	}

}



}	// mla

#endif
