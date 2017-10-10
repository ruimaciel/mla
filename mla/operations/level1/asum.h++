#ifndef MLA_OPERATIONS_LEVEL1_ASUM_HPP
#define MLA_OPERATIONS_LEVEL1_ASUM_HPP

#include <type_traits>
#include <cmath>

#include <mla/vector/all.h++>


namespace mla {

/**
 * Implements the level 1 BLAS functions
 * sum of absolute values
 * http://www.netlib.org/blas/#_level_1
 **/

template<typename Scalar, template<typename> class VectorStoragePolicyX>
Scalar 
asum(VectorStoragePolicyX<Scalar> &x)
{
	typename VectorStoragePolicyX<Scalar>::Cursor cursor = x.cursor();
	cursor.reset();

	Scalar absolute_sum = 0.0f;	// implicit value
	while(!cursor.at_end())
	{
		Scalar value = cursor.element();
		absolute_sum += abs(value);
		cursor.next();
	}

	return absolute_sum;
}


template<typename Scalar>
Scalar 
asum(mla::vector::Dense<Scalar> const &x)
{

	Scalar absolute_sum = 0.0f;	// implicit value
	for( auto it = x.data.begin(); it != x.data.end(); ++it)
	{
		Scalar value = *it;
		absolute_sum += abs(value);
	}

	return absolute_sum;
}


}	// mla

#endif
