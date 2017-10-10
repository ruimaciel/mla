#define BOOST_TEST_MODULE matrix

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/all.h++>
#include <mla/matrix/convert.h++>


typedef boost::mpl::list<
	mla::matrix::DenseRowMajor<float>,
	mla::matrix::DenseRowMajor<double>,
	mla::matrix::SparseDOK<float>,
	mla::matrix::SparseDOK<double>,
	mla::matrix::SparseCOO<float>,
	mla::matrix::SparseCOO<double>,
	mla::matrix::SparseCRS<float>,
	mla::matrix::SparseCRS<double>,
	mla::matrix::SparseCCS<float>,
	mla::matrix::SparseCCS<double>
> matrix_type_list;


BOOST_AUTO_TEST_SUITE(test_matrix)

BOOST_AUTO_TEST_CASE_TEMPLATE( matrix_convert, MatrixTypeTo, matrix_type_list )
{
	size_t matrix_size = 6;

	mla::matrix::DenseRowMajor<float> from(matrix_size, matrix_size);

	from(0,0) = (typename MatrixTypeTo::scalar_type)1;
	from(0,3) = (typename MatrixTypeTo::scalar_type)3;
	from(3,0) = (typename MatrixTypeTo::scalar_type)3;
	from(2,3) = (typename MatrixTypeTo::scalar_type)5;

	MatrixTypeTo to(matrix_size, matrix_size);


	mla::matrix::convert(from, to);

	BOOST_CHECK_EQUAL(from.rows(), to.rows());
	BOOST_CHECK_EQUAL(from.columns(), to.columns());

	for(size_t i = 0; i < matrix_size; i++)
	{
		for(size_t j = 0; i < matrix_size; i++)
		{
			BOOST_CHECK_CLOSE( from.getValue(i, j), to.getValue(i,j), 0.001f );
		}
	}
}


BOOST_AUTO_TEST_SUITE_END()

