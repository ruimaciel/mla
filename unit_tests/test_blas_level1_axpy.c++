#define BOOST_TEST_MODULE blas

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/vector/all.h++>

#include <mla/operations/level1/axpy.h++>


typedef boost::mpl::list<
	mla::vector::Dense<float>,
	mla::vector::Dense<double>,
	mla::vector::SparseCS<float>,
	mla::vector::SparseCS<double>
> vector_type_list;



BOOST_AUTO_TEST_SUITE(test_boost_level1)


BOOST_AUTO_TEST_CASE_TEMPLATE( blas_level1_test_axpy_zero, VectorType, vector_type_list )
{
	using namespace mla;

	VectorType x(3);
	vector::Dense<typename VectorType::scalar_type > y(3);

	y.setZero();
	
	typename VectorType::scalar_type value = 1.0f;
	axpy( value, x, y);

	for( unsigned int i = 0; i < y.size(); i++)
	{
		BOOST_CHECK_CLOSE(y.getValue(i), 0.0f, 0.001f);
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE( blas_level1_test_axpy_one, VectorType, vector_type_list )
{
	using namespace mla;

	VectorType x(3);
	x.setValue(0, 1.0f);
	x.setValue(1, 1.0f);
	x.setValue(2, 1.0f);

	vector::Dense<typename VectorType::scalar_type > y(3);
	y.setValue(0, 1.0f);
	y.setValue(1, 1.0f);
	y.setValue(2, 1.0f);

	typename VectorType::scalar_type value = 1.0f;
	axpy( value, x, y);

	for( unsigned int i = 0; i < y.size(); i++)
	{
		BOOST_CHECK_CLOSE(y.getValue(i), 2.0f, 0.001f);
	}
}


BOOST_AUTO_TEST_SUITE_END()

