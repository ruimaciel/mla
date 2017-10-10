#define BOOST_TEST_MODULE Test_VectorCursor

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/vector/all.h++>


typedef boost::mpl::list<
	mla::vector::Dense<float>,
	mla::vector::Dense<double>,
	mla::vector::SparseCS<float>,
	mla::vector::SparseCS<double>
> vector_type_list;



BOOST_AUTO_TEST_CASE_TEMPLATE( VectorCursor_iterate_through_dense, VectorType, vector_type_list )
{
	using namespace mla::vector;

	size_t size = 3;
	VectorType m(size);

	m.setZero();
	for(size_t i = 0; i < size; i++)
	{
		m.setValue(i, 1.0f*i);
	}

	typename VectorType::Cursor cursor(m);

	// checks if the cursor is able to return the right values
	for(size_t i = 0; i < size; i++)
	{
		BOOST_CHECK_CLOSE( cursor.element(), 1.0f*i, 1.0e-5 );
		cursor.next();
	}
}


