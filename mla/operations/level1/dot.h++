#ifndef MLA_OPERATIONS_LEVEL1_DOT_HPP
#define MLA_OPERATIONS_LEVEL1_DOT_HPP

#include <type_traits>

#include <mla/LAException.h++>

#include <mla/vector/all.h++>

namespace mla {

/**
 * Implements the level 1 BLAS functions
 * http://www.netlib.org/blas/#_reference_blas_version_3_5_0
 **/

template<typename Scalar, template<typename> class VectorStoragePolicyX, template<typename> class VectorStoragePolicyY>
Scalar 
dot(VectorStoragePolicyX<Scalar> &x, VectorStoragePolicyY<Scalar> &y)
{
	if( x.size() != y.size() )
	{
		throw LAException("level1::dot: incompatible vector size");
	}

	typename VectorStoragePolicyX<Scalar>::Cursor cursorx = x.cursor();
	typename VectorStoragePolicyY<Scalar>::Cursor cursory = y.cursor();

	cursorx.reset();
	cursory.reset();

	Scalar accumulator = 0.0f;

	size_t x_i = cursorx.current();
	size_t y_i = cursory.current();

	while( ! (cursorx.at_end() || cursory.at_end() ) )
	{
		if(x_i == y_i)
		{
			accumulator += cursorx.element()*cursory.element();
			cursorx.next();
			cursory.next();
		}
		else if (x_i > y_i)
		{
			cursory.next();
			y_i = cursory.current();
		}
		else
		{
			cursorx.next();
			x_i = cursorx.current();
		}
	}

	return accumulator;
}


template<typename Scalar, template<typename> class VectorStoragePolicyX>
Scalar 
dot(VectorStoragePolicyX<Scalar> &x, vector::Dense<Scalar> &y)
{
	if( x.size() != y.size() )
	{
		throw LAException("level1::dot: incompatible vector size");
	}

	typename VectorStoragePolicyX<Scalar>::Cursor cursorx = x.cursor();

	cursorx.reset();

	Scalar accumulator = 0.0f;


	while( !cursorx.at_end() )
	{
		size_t x_i = cursorx.current();
		accumulator += cursorx.element()*y.getValue(x_i);
		cursorx.next();
	}

	return accumulator;
}



}	// namespace mla

#endif
