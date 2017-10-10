#define BOOST_TEST_MODULE vector

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/vector/all.h++>
#include <mla/vector/convert.h++>


typedef boost::mpl::list<
	mla::vector::Dense<float>,
	mla::vector::Dense<double>,
	mla::vector::SparseCS<float>,
	mla::vector::SparseCS<double>
> vector_type_list;


using VectorTypeTo = mla::vector::SparseCS<float>;



BOOST_AUTO_TEST_SUITE(vector)


BOOST_AUTO_TEST_CASE_TEMPLATE( convert, VectorTypeFrom, vector_type_list )
{
	size_t vector_size = 6;

	VectorTypeFrom from(vector_size);
	
	from.setValue(1, 1.0f);
	from.setValue(3, 3.0f);
	from.setValue(5, 5.0f);

	//VectorTypeTo to(vector_size);	// g++ throws segfault
	mla::vector::SparseCS<float> to(vector_size);


	mla::vector::convert(from, to, 0.001f);

	BOOST_CHECK_EQUAL(to.data.values.size(), 3);
	BOOST_CHECK_EQUAL(to.data.values[0], 1.0f);
	BOOST_CHECK_EQUAL(to.data.values[1], 3.0f);
	BOOST_CHECK_EQUAL(to.data.values[2], 5.0f);
	BOOST_CHECK_EQUAL(to.data.column_index[0], 1);
	BOOST_CHECK_EQUAL(to.data.column_index[1], 3);
	BOOST_CHECK_EQUAL(to.data.column_index[2], 5);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( convert_ignore_nulls, VectorTypeFrom, vector_type_list )
{
	/**
	 * Test intended to check if convert actually ignores values lower than the 
	 * established non-null limit.
	 */

	size_t vector_size = 6;
	float non_null_limit = 0.001f;	// the non-null limit for this case

	VectorTypeFrom from(vector_size);
	
	from.setValue(1, 1.0f);
	from.setValue(2, non_null_limit/2);
	from.setValue(3, 3.0f);
	from.setValue(5, 5.0f);

	//VectorTypeTo to(vector_size);	// g++ throws segfault
	mla::vector::SparseCS<float> to(vector_size);


	mla::vector::convert(from, to, non_null_limit);

	BOOST_CHECK_EQUAL(to.data.values.size(), 3);
	BOOST_CHECK_EQUAL(to.data.values[0], 1.0f);
	BOOST_CHECK_EQUAL(to.data.values[1], 3.0f);
	BOOST_CHECK_EQUAL(to.data.values[2], 5.0f);
	BOOST_CHECK_EQUAL(to.data.column_index[0], 1);
	BOOST_CHECK_EQUAL(to.data.column_index[1], 3);
	BOOST_CHECK_EQUAL(to.data.column_index[2], 5);
}



BOOST_AUTO_TEST_SUITE_END()

