#define BOOST_TEST_MODULE matrix

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/all.h++>


typedef boost::mpl::list<
	mla::matrix::DenseRowMajor<float>,
	mla::matrix::DenseRowMajor<double>,
	mla::matrix::Diagonal<float>,
	mla::matrix::Diagonal<double>,
	mla::matrix::SparseDOK<float>,
	mla::matrix::SparseDOK<double>,
	mla::matrix::SparseCRS<float>,
	mla::matrix::SparseCRS<double>
> matrix_type_list;


BOOST_AUTO_TEST_SUITE(test_matrix)


BOOST_AUTO_TEST_CASE_TEMPLATE( Declaration, MatrixType, matrix_type_list )
{
	size_t row_size = 3;
	size_t column_size = 3;
	
	MatrixType m(row_size,column_size);

	BOOST_CHECK_EQUAL(m.rows(), row_size);
	BOOST_CHECK_EQUAL(m.columns(), column_size);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( Clearing, MatrixType, matrix_type_list )
{
	size_t row_size = 3;
	size_t column_size = 3;

	MatrixType m(row_size,column_size);

	m.setZero();

	for(unsigned int i = 0; i < row_size; i++)
	{
		for(unsigned int j = 0; j < column_size; j++)
		{
			BOOST_CHECK_EQUAL(m.getValue(i, j), (typename MatrixType::scalar_type)0);
		}
	}
	
}


BOOST_AUTO_TEST_CASE_TEMPLATE( set_eye_square, MatrixType, matrix_type_list )
{
	size_t row_size = 3;
	size_t column_size = 3;

	MatrixType m(row_size,column_size);

	m.setEye();

	for(unsigned int i = 0; i < row_size; i++)
	{
		for(unsigned int j = 0; j < column_size; j++)
		{
			BOOST_CHECK_EQUAL(m.getValue(i, j), (typename MatrixType::scalar_type)(i==j?1:0) );
		}
	}
	
}



BOOST_AUTO_TEST_CASE_TEMPLATE( set_eye_rectangular_row, MatrixType, matrix_type_list )
{
	size_t row_size = 3;
	size_t column_size = 5;

	MatrixType m(row_size,column_size);

	m.setEye();

	for(unsigned int i = 0; i < row_size; i++)
	{
		for(unsigned int j = 0; j < column_size; j++)
		{
			BOOST_CHECK_EQUAL(m.getValue(i, j), (typename MatrixType::scalar_type)(i==j?1:0) );
		}
	}
	
}


BOOST_AUTO_TEST_CASE_TEMPLATE( set_eye_rectangular_column, MatrixType, matrix_type_list )
{
	size_t row_size = 5;
	size_t column_size = 3;

	MatrixType m(row_size, column_size);

	m.setEye();

	for(unsigned int i = 0; i < row_size; i++)
	{
		for(unsigned int j = 0; j < column_size; j++)
		{
			BOOST_CHECK_EQUAL(m.getValue(i, j), (typename MatrixType::scalar_type)(i==j?1:0) );
		}
	}
	
}



BOOST_AUTO_TEST_CASE_TEMPLATE( reference_assigning_one, MatrixType, matrix_type_list )
{
	if( mla::matrix::Traits<MatrixType>::is_writeable() )
	{
		MatrixType m(3,3);
		m.setZero();

		m(1,1) = 1.0;	

		BOOST_CHECK_EQUAL(m.getValue(0, 0), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(0, 1), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(0, 2), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(1, 0), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(1, 1), 1.0);
		BOOST_CHECK_EQUAL(m.getValue(1, 2), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(2, 0), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(2, 1), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(2, 2), 0.0);
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE( reference_assigning_same_row_ordered, MatrixType, matrix_type_list )
{
	if( mla::matrix::Traits<MatrixType>::is_writeable() )
	{
		MatrixType m(3,3);
		m.setZero();

		m(1,0) = 1.0;
		m(1,1) = 1.0;

		BOOST_CHECK_EQUAL(m.getValue(0, 0), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(0, 1), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(0, 2), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(1, 0), 1.0);
		BOOST_CHECK_EQUAL(m.getValue(1, 1), 1.0);
		BOOST_CHECK_EQUAL(m.getValue(1, 2), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(2, 0), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(2, 1), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(2, 2), 0.0);
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE( reference_assigning_same_row_unordered, MatrixType, matrix_type_list )
{
	if( mla::matrix::Traits<MatrixType>::is_writeable() )
	{
		MatrixType m(3,3);
		m.setZero();

		m(1,1) = 1.0;
		m(1,0) = 1.0;

		BOOST_CHECK_EQUAL(m.getValue(0, 0), 0);
		BOOST_CHECK_EQUAL(m.getValue(0, 1), 0);
		BOOST_CHECK_EQUAL(m.getValue(0, 2), 0);
		BOOST_CHECK_EQUAL(m.getValue(1, 0), 1);
		BOOST_CHECK_EQUAL(m.getValue(1, 1), 1);
		BOOST_CHECK_EQUAL(m.getValue(1, 2), 0);
		BOOST_CHECK_EQUAL(m.getValue(2, 0), 0);
		BOOST_CHECK_EQUAL(m.getValue(2, 1), 0);
		BOOST_CHECK_EQUAL(m.getValue(2, 2), 0);
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE( reference_assigning_different_row_ordered, MatrixType, matrix_type_list )
{
	if( mla::matrix::Traits<MatrixType>::is_writeable() )
	{
		MatrixType m(3,3);
		m.setZero();

		m(1,1) = 1.0;
		m(2,1) = 1.0;

		BOOST_CHECK_EQUAL(m.getValue(0, 0), 0);
		BOOST_CHECK_EQUAL(m.getValue(0, 1), 0);
		BOOST_CHECK_EQUAL(m.getValue(0, 2), 0);
		BOOST_CHECK_EQUAL(m.getValue(1, 0), 0);
		BOOST_CHECK_EQUAL(m.getValue(1, 1), 1);
		BOOST_CHECK_EQUAL(m.getValue(1, 2), 0);
		BOOST_CHECK_EQUAL(m.getValue(2, 0), 0);
		BOOST_CHECK_EQUAL(m.getValue(2, 1), 1);
		BOOST_CHECK_EQUAL(m.getValue(2, 2), 0);
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE( value_setting, MatrixType, matrix_type_list )
{
	if( mla::matrix::Traits<MatrixType>::is_writeable() )
	{
		MatrixType m(3,3);
		m.setZero();

		m.setValue( 1, 1, (typename MatrixType::scalar_type)1);	

		BOOST_CHECK_EQUAL(m.getValue(1, 1), 1.0);
		BOOST_CHECK_EQUAL(m.getValue(0, 0), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(0, 1), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(0, 2), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(1, 0), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(1, 1), 1.0);
		BOOST_CHECK_EQUAL(m.getValue(1, 2), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(2, 0), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(2, 1), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(2, 2), 0.0);
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE( matrix_assigning, MatrixType, matrix_type_list )
{
	if( mla::matrix::Traits<MatrixType>::is_writeable() )
	{
		MatrixType m(3,3), n;
		m.setZero();

		m.setValue( 1, 1, 1.0);	

		n = m;

		BOOST_CHECK_EQUAL(n.getValue(1, 1), 1.0);
		BOOST_CHECK_EQUAL(n.getValue(0, 0), 0.0);
		BOOST_CHECK_EQUAL(n.getValue(0, 1), 0.0);
		BOOST_CHECK_EQUAL(n.getValue(0, 2), 0.0);
		BOOST_CHECK_EQUAL(n.getValue(1, 0), 0.0);
		BOOST_CHECK_EQUAL(n.getValue(1, 1), 1.0);
		BOOST_CHECK_EQUAL(n.getValue(1, 2), 0.0);
		BOOST_CHECK_EQUAL(n.getValue(2, 0), 0.0);
		BOOST_CHECK_EQUAL(n.getValue(2, 1), 0.0);
		BOOST_CHECK_EQUAL(m.getValue(2, 2), 0.0);
	}
}


BOOST_AUTO_TEST_SUITE_END()


