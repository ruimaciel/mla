#define BOOST_TEST_MODULE blas

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/all.h++>
#include <mla/vector/all.h++>

#include <mla/operations/level2/gemv.h++>


using Scalar = float;

typedef boost::mpl::list<
	mla::matrix::DenseRowMajor<Scalar>,
	mla::matrix::Diagonal<Scalar>,
	mla::matrix::SparseDOK<Scalar>,
	mla::matrix::SparseCRS<Scalar>
> matrix_type_list;


using MatrixType = mla::matrix::DenseRowMajor<Scalar>;


typedef boost::mpl::list<
	mla::vector::Dense<Scalar>,
	mla::vector::SparseCS<Scalar>
> vector_type_list;



BOOST_AUTO_TEST_SUITE(test_boost_level2)


BOOST_AUTO_TEST_CASE_TEMPLATE( blas_level2_test_gemv_zero, VectorType, vector_type_list )
{
	using namespace mla;

	Scalar alpha = 1.0f;

	//MatrixType A(3,3);
	size_t length = 3;
	mla::matrix::DenseRowMajor<Scalar> A( length, length);

	A(0,0) = 1.0f;
	A(1,1) = 1.0f;
	A(2,2) = 1.0f;

	VectorType x(length);
	x.setValue(0, 1.0f);
	x.setValue(1, 1.0f);
	x.setValue(2, 1.0f);

	Scalar beta = 1.0f;

	VectorType y(length);
	y.setValue(0, 1.0f);
	y.setValue(1, 1.0f);
	y.setValue(2, 1.0f);

	gemv(alpha, A, x, beta, y);


	for(size_t i = 0; i < length; i++)
	{
		BOOST_CHECK_CLOSE( y.getValue(i), 2.0f, 1.0e-5 );
	}
}



BOOST_AUTO_TEST_SUITE_END()

