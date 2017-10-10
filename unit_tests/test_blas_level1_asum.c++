#define BOOST_TEST_MODULE blas

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/vector/all.h++>

#include <mla/operations/level1/asum.h++>


typedef boost::mpl::list<
	mla::vector::Dense<float>,
	mla::vector::Dense<double>,
	mla::vector::SparseCS<float>,
	mla::vector::SparseCS<double>
> vector_type_list;



BOOST_AUTO_TEST_SUITE(test_boost_level1)


BOOST_AUTO_TEST_CASE_TEMPLATE( blas_level1_test_asum_one_dense, VectorType, vector_type_list )
{
	using namespace mla;

	VectorType x(3);
	x.setValue(0, 1.0f);
	x.setValue(1, 1.0f);
	x.setValue(2, 1.0f);

	typename VectorType::scalar_type value;
	value = asum(x);

	BOOST_CHECK_CLOSE( value, 3.0f, 1.0e-5 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( blas_level1_test_asum_one_sparse, VectorType, vector_type_list )
{
	using namespace mla;

	VectorType x(5);
	x.setValue(0, 1.0f);
	x.setValue(2, 1.0f);
	x.setValue(4, 1.0f);

	typename VectorType::scalar_type value;
	value = asum(x);

	BOOST_CHECK_CLOSE( value, 3.0f, 1.0e-5 );
}


BOOST_AUTO_TEST_SUITE_END()

