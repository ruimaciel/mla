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


using VectorTypeFrom = mla::vector::Dense<float>;



BOOST_AUTO_TEST_SUITE(vector)


BOOST_AUTO_TEST_CASE_TEMPLATE( convert, VectorTypeTo, vector_type_list )
{
	size_t vector_size = 6;

	//VectorTypeFrom from(vector_size);	// g++ throws segfault
	mla::vector::Dense<float> from(vector_size);

	from.setValue(1, 1.0f);
	from.setValue(3, 3.0f);
	from.setValue(5, 5.0f);

	VectorTypeTo to(vector_size);


	mla::vector::convert(from, to);

	BOOST_CHECK_EQUAL(from.size(), to.size());

	for(size_t i = 0; i < vector_size; i++)
	{
		BOOST_CHECK_EQUAL(from.getValue(i), to.getValue(i));
	}
}



BOOST_AUTO_TEST_SUITE_END()

