#ifndef MLA_MATRIX_CONVERT_HPP
#define MLA_MATRIX_CONVERT_HPP


#include <cmath>        // std::abs

#include <mla/LAException.h++>

#include <mla/matrix/all.h++>

namespace mla
{
namespace matrix
{


/**
 * Generic routine to convert between any Matrix class, using only the generic interface
 *@param from	the origin matrix, which is to be converted to another format
 *@param to	the destination matrix, which is to be converted to another format
 *@param interpret_as_zero_limit	all elements that are below this number will be interpreted as null entries
 */
template<typename FromScalar, template <typename> class FromMatrix, typename ToScalar, template <typename> class ToMatrix>
void
convert(FromMatrix<FromScalar> &from, ToMatrix<ToScalar> &to, double interpret_as_zero_limit = 0.0f)
{
	auto from_cursor = from.cursor();

	interpret_as_zero_limit = std::abs(interpret_as_zero_limit);

	// set the matrix size
	to.resize( from.rows(), from.columns() );

	// iterate over each row
	size_t i;
	size_t j;
	while( !from_cursor.at_end_of_rows() )
	{
		i = from_cursor.current_row();
		
		// iterate over each element in the row
		while( !from_cursor.at_end_of_current_row() )
		{
			ToScalar value = (ToScalar)from_cursor.element();

			if( value > interpret_as_zero_limit)
			{
				j = from_cursor.current_column();

				to(i,j) = (ToScalar)from_cursor.element();
			}

			from_cursor.increment_column();
		}

		// move to next non-null element in the next row
		from_cursor.start_next_row_nn();
	}
}


/**
 * Partial template specialization from a SparseDOK matrix to any other matrix
 *@param from	the origin matrix, which is a SparseDOK data structure
 *@param to	the destination matrix, which is to be converted to another format
 *@param interpret_as_zero_limit	all elements that are below this number will be interpreted as null entries
 */
template<typename FromScalar, typename ToScalar, template <typename> class ToMatrix>
void
convert(SparseDOK<FromScalar> &from, ToMatrix<ToScalar> &to, double interpret_as_zero_limit = 0.0f)
{
	interpret_as_zero_limit = std::abs(interpret_as_zero_limit);
	
	size_t row, column;
	for(typename std::map< size_t, FromScalar>::const_iterator iter = from.data.key_value_map.begin(); iter != from.data.key_value_map.end(); iter++)
	{
		const ToScalar &value =  iter->second;
		if( std::abs(value) > interpret_as_zero_limit )
		{
			row = iter->first % from.columns();
			column = iter->first - row*from.columns();
			to(row, column) = (ToScalar)value;
		}
	}
}


/**
 * Partial template specialization from a SparseCOO matrix to any other matrix
 *@param from	the origin matrix, which is a SparseCOO data structure
 *@param to	the destination matrix, which is to be converted to another format
 *@param interpret_as_zero_limit	all elements that are below this number will be interpreted as null entries
 */
template<typename FromScalar, typename ToScalar, template <typename> class ToMatrix>
void
convert(SparseCOO<FromScalar> &from, ToMatrix<ToScalar> &to, double interpret_as_zero_limit = 0.0f)
{
	interpret_as_zero_limit = std::abs(interpret_as_zero_limit);
	
	size_t row, column;
	for(typename std::list< typename SparseCOO<FromScalar>::Data::Coordinate>::const_iterator iter = from.data.coordinate_list.begin(); iter != from.data.coordinate_list.end(); iter++)
	{
		ToScalar const value =  iter->value;
		if( std::abs(value) > interpret_as_zero_limit )
		{
			row = iter->i;
			column = iter->j;
			to(row, column) = (ToScalar)value;
		}
	}
}


}	// namespace mla::matrix

}	// namespace mla

#endif

