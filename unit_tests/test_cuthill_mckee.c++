#define BOOST_TEST_MODULE matrix

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/all.h++>
#include <mla/matrix/convert.h++>
#include <mla/operations/cuthill_mckee.h++>


typedef boost::mpl::list<
//	mla::matrix::DenseRowMajor<float>,
//	mla::matrix::DenseRowMajor<double>,
//	mla::matrix::SparseCRS<float>,
//	mla::matrix::SparseCRS<double>,
	mla::matrix::SparseDOK<float>,
	mla::matrix::SparseDOK<double>
> matrix_type_list;


BOOST_AUTO_TEST_SUITE(test_operations)


BOOST_AUTO_TEST_CASE_TEMPLATE( matrix_cuthill_mckee_arrowhead, MatrixType, matrix_type_list )
{
	size_t matrix_size = 6;

	typedef typename MatrixType::scalar_type Scalar;

	MatrixType A(matrix_size, matrix_size);

	// sets an arrowhead matrix
	for(unsigned int i = 0; i < matrix_size; i++)
	{
		A.setValue( i, i, (Scalar)1.0f );
		A.setValue( 0, i, (Scalar)1.0f );
		A.setValue( i, 0, (Scalar)1.0f );
	}

	std::vector<size_t> indices = mla::cuthill_mckee(A);


	// the indices should be the same size as the matrix
	BOOST_CHECK_EQUAL( indices.size(), matrix_size );

	// Cuthill-McKee reorders matrices as arrowhead matrices
	for(size_t i = 0; i < matrix_size; i++)
	{
		BOOST_CHECK_EQUAL( indices[i], i );
	}

}


BOOST_AUTO_TEST_SUITE_END()

