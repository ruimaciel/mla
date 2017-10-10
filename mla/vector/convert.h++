#ifndef MLA_VECTOR_CONVERT_HPP
#define MLA_VECTOR_CONVERT_HPP


#include <cmath>        // std::abs

#include <mla/LAException.h++>

#include <mla/vector/all.h++>

namespace mla
{
namespace vector
{


template<typename FromScalar, template <typename> class FromVector, typename ToScalar, template <typename> class ToVector>
void
convert(FromVector<FromScalar> &from, ToVector<ToScalar> &to, double interpret_as_zero_limit = 0.0f)
{
	auto from_cursor = from.cursor();

	interpret_as_zero_limit = std::abs(interpret_as_zero_limit);

	// set the vector size
	to.resize( from.size() );

	// iterate over each row
	size_t i;
	while( !from_cursor.at_end() )
	{
		// iterate over each element in the row
		ToScalar value = (ToScalar)from_cursor.element();
		if( value > interpret_as_zero_limit)
		{
			i = from_cursor.current();
			to.setValue(i, value);
		}
			
		from_cursor.next();
	}
}


}	// namespace mla::vector

}	// namespace mla

#endif
