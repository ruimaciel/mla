#ifndef MLA_OPERATIONS_LEVEL1_SCALE_HPP
#define MLA_OPERATIONS_LEVEL1_SCALE_HPP

#include <type_traits>

#include <mla/vector/all.h++>


namespace mla {

/**
 * Implements the level 1 BLAS functions
 * {x} := a*{x}
 * http://www.netlib.org/blas/#_reference_blas_version_3_5_0
 **/

template<typename Scalar, template<typename> class VectorStoragePolicyX>
void 
scale(Scalar a, VectorStoragePolicyX<Scalar> &x)
{
	typename VectorStoragePolicyX<Scalar>::Cursor cursor = x.cursor();
	cursor.reset();

	while(!cursor.at_end())
	{
		size_t i = cursor.current();
		x.setValue(i, cursor.element()*a);
		cursor.next();
	}
}


}	// mla

#endif
