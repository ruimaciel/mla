#define BOOST_TEST_MODULE matrix

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/all.h++>
#include <mla/vector/all.h++>
#include <mla/matrix/convert.h++>

#include <mla/solvers/CG.h++>


typedef boost::mpl::list<
	mla::matrix::DenseRowMajor<float>,
	mla::matrix::DenseRowMajor<double>,
	mla::matrix::SparseCRS<float>,
	mla::matrix::SparseCRS<double>
> matrix_type_list;


BOOST_AUTO_TEST_SUITE(test_solvers)


BOOST_AUTO_TEST_CASE_TEMPLATE( solve_unit_1, MatrixType, matrix_type_list )
{
	size_t matrix_size = 6;

	typedef typename MatrixType::scalar_type Scalar;

	MatrixType A(matrix_size, matrix_size);
	A.setEye();

	mla::vector::Dense<Scalar> x(matrix_size), b(matrix_size);

	Scalar value = 1.0f;

	for(unsigned int i = 0; i < matrix_size; i++)
	{
		b.setValue( i, value);
	}



	mla::cg(A, x, b, 0.01f, 6);

	for(unsigned int i = 0; i < matrix_size; i++)
	{
		BOOST_CHECK_CLOSE( x.getValue(i), 1.0f/b.getValue(i), 0.001f);
	}
}



BOOST_AUTO_TEST_CASE_TEMPLATE( solve_unit_2, MatrixType, matrix_type_list )
{
	size_t matrix_size = 6;

	typedef typename MatrixType::scalar_type Scalar;

	MatrixType A(matrix_size, matrix_size);
	A.setEye();

	mla::vector::Dense<Scalar> x(matrix_size), b(matrix_size);

	Scalar value = 2.0f;

	for(unsigned int i = 0; i < matrix_size; i++)
	{
		b.setValue( i, value);
	}


	mla::cg(A, x, b, 0.01f, 6);

	for(unsigned int i = 0; i < matrix_size; i++)
	{
		BOOST_CHECK_CLOSE( x.getValue(i), b.getValue(i), 0.001f);
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE( solve_unit_3, MatrixType, matrix_type_list )
{
	size_t matrix_size = 6;

	typedef typename MatrixType::scalar_type Scalar;

	mla::matrix::DenseRowMajor<Scalar> from(matrix_size, matrix_size);
	mla::vector::Dense<Scalar> x(matrix_size), b(matrix_size);

	Scalar value = 2.0f;

	for(unsigned int i = 0; i < matrix_size; i++)
	{
		from.setValue( i, i, value);
		b.setValue( i, (Scalar)1.0f);
	}

	MatrixType A(matrix_size, matrix_size);
	mla::matrix::convert(from, A);

	mla::cg(A, x, b, 0.01f, 6);

	for(unsigned int i = 0; i < matrix_size; i++)
	{
		BOOST_CHECK_CLOSE( x.getValue(i), 1.0f/value, 0.001f);
	}
}



BOOST_AUTO_TEST_SUITE_END()

