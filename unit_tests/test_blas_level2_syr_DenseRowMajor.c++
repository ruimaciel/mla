#define BOOST_TEST_MODULE blas

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/all.h++>
#include <mla/vector/all.h++>

#include <mla/operations/level2/syr.h++>


using Scalar = float;

typedef boost::mpl::list<
	mla::matrix::DenseRowMajor<Scalar>,
	mla::matrix::Diagonal<Scalar>,
	mla::matrix::SparseDOK<Scalar>,
	mla::matrix::SparseCRS<Scalar>
> matrix_type_list;


template<typename T> using MatrixType = mla::matrix::DenseRowMajor<T>;


typedef boost::mpl::list<
	mla::vector::Dense<Scalar>,
	mla::vector::SparseCS<Scalar>
> vector_type_list;



BOOST_AUTO_TEST_SUITE(test_boost_level2)


BOOST_AUTO_TEST_CASE_TEMPLATE( blas_level2_test_syr_DenseRowMajor_zero, VectorType, vector_type_list )
{
	using namespace mla;

	Scalar alpha = 1.0f;

	size_t length = 3;
	MatrixType<Scalar> A( length, length);

	VectorType x(length);
	x.setValue(0, 1.0f);
	x.setValue(1, 1.0f);
	x.setValue(2, 1.0f);


	syr(alpha, x, A);


	for(size_t i = 0; i < length; i++)
	{
		for(size_t j = 0; j < length; j++)
		{
			BOOST_CHECK_CLOSE( A.getValue(i,j), 1.0f, 1.0e-5 );
		}
	}
}



BOOST_AUTO_TEST_SUITE_END()

