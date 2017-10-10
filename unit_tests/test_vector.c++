#define BOOST_TEST_MODULE vector

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



BOOST_AUTO_TEST_SUITE(vector)


BOOST_AUTO_TEST_CASE_TEMPLATE( Declaration, VectorType, vector_type_list )
{
	VectorType v(3);

	BOOST_CHECK_EQUAL(v.size(), 3);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( Clearing, VectorType, vector_type_list )
{
	VectorType v(3);
	v.setZero();

	for(unsigned int i = 0; i < 3; i++)
	{
		BOOST_CHECK_EQUAL(v.getValue(i), 0.0);
	}
	
}


BOOST_AUTO_TEST_CASE_TEMPLATE( Assigning, VectorType, vector_type_list )
{
	VectorType v(3);
	v.setZero();

	v.setValue(1, 1.0);	

	BOOST_CHECK_EQUAL(v.getValue(0), 0.0);
	BOOST_CHECK_EQUAL(v.getValue(1), 1.0);
	BOOST_CHECK_EQUAL(v.getValue(2), 0.0);
}


BOOST_AUTO_TEST_SUITE_END()

