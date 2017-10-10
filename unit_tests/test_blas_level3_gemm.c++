#define BOOST_TEST_MODULE matrix

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/all.h++>
#include <mla/vector/all.h++>
#include <mla/matrix/convert.h++>

#include <mla/operations/level3/gemm.h++>


typedef boost::mpl::list<
	mla::matrix::DenseRowMajor<float>,
	mla::matrix::DenseRowMajor<double>
> matrix_type_list;


BOOST_AUTO_TEST_SUITE(test_operations)


BOOST_AUTO_TEST_CASE_TEMPLATE( matrix_matrix_multiply_identity, MatrixType, matrix_type_list )
{
	size_t matrix_size = 6;

	typedef typename MatrixType::scalar_type Scalar;

	mla::matrix::DenseRowMajor<Scalar> from(matrix_size, matrix_size);

	for(unsigned int i = 0; i < matrix_size; i++)
	{
		from.setValue( i, i, (Scalar)1.0f );
	}

	MatrixType A(matrix_size, matrix_size);
	MatrixType B(matrix_size, matrix_size);
	MatrixType C(matrix_size, matrix_size);

	mla::matrix::convert(from, A);
	mla::matrix::convert(from, B);
	C.setZero();

	Scalar const alfa = 1.0;
	Scalar const beta = 1.0;
	mla::gemm(alfa, A, B, beta, C);

	for(unsigned int i = 0; i < matrix_size; i++)
	{
		for(unsigned int j = 0; j < matrix_size; j++)
		{
			BOOST_CHECK_CLOSE( C.getValue(i,j), A.getValue(i,j), 0.001f);
		}
	}
}


BOOST_AUTO_TEST_SUITE_END()

