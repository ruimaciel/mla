#define BOOST_TEST_MODULE matrix

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/all.h++>
#include <mla/vector/all.h++>
#include <mla/matrix/convert.h++>

#include <mla/operations/level2/gemv.h++>


typedef boost::mpl::list<
	mla::matrix::DenseRowMajor<float>,
	mla::matrix::DenseRowMajor<double>,
	mla::matrix::SparseCRS<float>,
	mla::matrix::SparseCRS<double>
> matrix_type_list;


BOOST_AUTO_TEST_SUITE(test_operations)


BOOST_AUTO_TEST_CASE_TEMPLATE( matrix_multiply, MatrixType, matrix_type_list )
{
	size_t matrix_size = 6;

	typedef typename MatrixType::scalar_type Scalar;

	MatrixType A(matrix_size, matrix_size);
	A.setEye();

	mla::vector::Dense<Scalar> x(matrix_size), y(matrix_size);

	for(unsigned int i = 0; i < matrix_size; i++)
	{
		x.setValue( i, (Scalar)1.0f);
	}


	Scalar const alfa = 1.0;
	Scalar const beta = 1.0;
	mla::gemv(alfa, A, x, beta, y);

	for(unsigned int i = 0; i < matrix_size; i++)
	{
		BOOST_CHECK_CLOSE( x.getValue(i), y.getValue(i), 0.001f);
	}
}


BOOST_AUTO_TEST_SUITE_END()

