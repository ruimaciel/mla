#ifndef MLA_OPERATIONS_LEVEL1_AXPY_HPP
#define MLA_OPERATIONS_LEVEL1_AXPY_HPP

#include <type_traits>

#include <mla/vector/all.h++>

namespace mla {

/**
 * Implements the level 1 BLAS functions
 * {y} := a*{x} + {y}
 *
 * http://www.netlib.org/blas/#_reference_blas_version_3_5_0
 * https://software.intel.com/en-us/node/468394
 **/

template<typename Scalar, template<typename> class VectorStoragePolicyX, template<typename> class VectorStoragePolicyY>
void 
axpy(Scalar const a, VectorStoragePolicyX<Scalar> &x, VectorStoragePolicyY<Scalar> &y)
{
	typename VectorStoragePolicyX<Scalar>::Cursor cursor = x.cursor();
	cursor.reset();

	while(!cursor.at_end())
	{
		size_t i = cursor.current();
		Scalar temp = y.getValue(i);
		y.setValue(i, temp + cursor.element()*a);
		cursor.next();
	}
}



// explicit instantiation

template<typename Scalar>
void 
axpy(Scalar const a, mla::vector::Dense<Scalar> &x, mla::vector::Dense<Scalar> &y)
{
	for(size_t i = 0; i < y.size(); i++)
	{
		y[i] += a*x[i];
	}
}


}	// mla

#endif
