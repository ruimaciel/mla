#define BOOST_TEST_MODULE blas

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/vector/all.h++>

#include <mla/operations/level1/scale.h++>


typedef boost::mpl::list<
	mla::vector::Dense<float>,
	mla::vector::Dense<double>,
	mla::vector::SparseCS<float>,
	mla::vector::SparseCS<double>
> vector_type_list;



BOOST_AUTO_TEST_SUITE(test_boost_level1)


BOOST_AUTO_TEST_CASE_TEMPLATE( blas_level1_test_scale_one, VectorType, vector_type_list )
{
	using namespace mla;

	VectorType x(3);
	x.setValue(0, 1.0f);
	x.setValue(1, 1.0f);
	x.setValue(2, 1.0f);

	typename VectorType::scalar_type value = 1.0f;
	scale( value, x);

	for( unsigned int i = 0; i < x.size(); i++)
	{
		BOOST_CHECK_CLOSE(x.getValue(i), 1.0f, 1e-6f);
	}
}



BOOST_AUTO_TEST_SUITE_END()

