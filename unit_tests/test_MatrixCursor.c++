#define BOOST_TEST_MODULE MatrixCursor

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/DenseRowMajor.h++>
#include <mla/matrix/Diagonal.h++>
#include <mla/matrix/SparseDOK.h++>
#include <mla/matrix/SparseCRS.h++>


typedef boost::mpl::list<
	mla::matrix::DenseRowMajor<float>,
	mla::matrix::DenseRowMajor<double>,
	mla::matrix::Diagonal<float>,
	mla::matrix::Diagonal<double>,
	mla::matrix::SparseCRS<float>,
	mla::matrix::SparseCRS<double>
> matrix_type_list;



BOOST_AUTO_TEST_SUITE(test_MatrixCursor)


BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixCursor_increment_row, MatrixType, matrix_type_list )
{
	using namespace mla::matrix;

	MatrixType m(3,3);
	m.setZero();

	typename MatrixType::Cursor cursor(m);

	// tests row incrementing
	BOOST_CHECK_EQUAL( cursor.current_row(), 0 );
	cursor.increment_row();
	BOOST_CHECK_EQUAL( cursor.current_row(), 1 );

}


BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixCursor_increment_column, MatrixType, matrix_type_list )
{
	using namespace mla::matrix;

	MatrixType m(3,3);
	m.setZero();

	typename MatrixType::Cursor cursor(m);

	// tests column incrementing
	BOOST_CHECK_EQUAL( cursor.current_column(), 0 );
	cursor.increment_column();
	BOOST_CHECK_EQUAL( cursor.current_column(), 1 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixCursor_reset, MatrixType, matrix_type_list )
{
	using namespace mla::matrix;

	MatrixType m(3,3);
	m.setZero();

	typename MatrixType::Cursor cursor(m);

	// tests row incrementing, and reseting the position
	cursor.increment_row();
	cursor.reset();
	BOOST_CHECK_EQUAL( cursor.current_row(), 0 );

	// tests column incrementing, and reseting the position
	cursor.increment_column();
	cursor.reset();
	BOOST_CHECK_EQUAL( cursor.current_column(), 0 );

	// tests both
	cursor.increment_row();
	cursor.increment_column();
	cursor.reset();
	BOOST_CHECK_EQUAL( cursor.current_row(), 0 );
	BOOST_CHECK_EQUAL( cursor.current_column(), 0 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixCursor_check_start, MatrixType, matrix_type_list )
{
	using namespace mla::matrix;

	MatrixType m(3,3);
	m.setZero();

	typename MatrixType::Cursor cursor(m);

	BOOST_CHECK( cursor.at_beginning_of_column() );
	BOOST_CHECK( cursor.at_beginning_of_row() );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixCursor_check_end_of_column, MatrixType, matrix_type_list )
{
	using namespace mla::matrix;

	MatrixType m(3,3);
	m.setZero();

	typename MatrixType::Cursor cursor(m);

	cursor.increment_column();
	cursor.increment_column();
	BOOST_CHECK( !cursor.at_end_of_columns() );
	cursor.increment_column();
	BOOST_CHECK( cursor.at_end_of_columns() );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixCursor_check_end_of_row, MatrixType, matrix_type_list )
{
	using namespace mla::matrix;

	MatrixType m(3,3);
	m.setZero();

	typename MatrixType::Cursor cursor(m);

	cursor.increment_row();
	cursor.increment_row();
	BOOST_CHECK( !cursor.at_end_of_rows() );
	cursor.increment_row();
	BOOST_CHECK( cursor.at_end_of_rows() );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixCursor_iterator_moving_test, MatrixType, matrix_type_list )
{
	using namespace mla::matrix;

	MatrixType m(3,3);
	m.setZero();
	m(0,0) = 1;
	m(1,1) = 2;
	m(2,2) = 3;

	typename MatrixType::Cursor cursor(m);

	cursor.reset();
	BOOST_CHECK_EQUAL( cursor.element(), 1.0 );
	cursor.increment_row();
	BOOST_CHECK_EQUAL( cursor.element(), 0.0 );
	cursor.increment_column();
	BOOST_CHECK_EQUAL( cursor.element(), 2.0 );
	cursor.increment_row();
	BOOST_CHECK_EQUAL( cursor.element(), 0.0 );
	cursor.increment_column();
	BOOST_CHECK_EQUAL( cursor.element(), 3.0 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixCursor_iterator_test_end_check, MatrixType, matrix_type_list )
{
	using namespace mla::matrix;

	MatrixType m(3,3);
	m.setZero();
	m(0,0) = 1;
	m(1,1) = 2;
	m(2,2) = 3;

	typename MatrixType::Cursor cursor(m);

	cursor.reset();
	for(size_t n = 0; n < 2; n++)
	{
		cursor.increment_row();
		cursor.increment_column();
		BOOST_CHECK_MESSAGE( !cursor.at_end_of_current_column(), "column(" << cursor.current_column() << ") is the end of a (" << m.rows() << "," << m.columns() << ") matrix" );
		BOOST_CHECK_MESSAGE( !cursor.at_end_of_current_row(), "row(" << cursor.current_row() << ") is the end of a (" << m.rows() << "," << m.columns() << ") matrix" );
	}

	// move past the last row and column
	cursor.increment_row();
	cursor.increment_column();
	BOOST_CHECK_MESSAGE( cursor.at_end_of_current_column(), "column(" << cursor.current_column() << ") is not the end of a (" << m.rows() << "," << m.columns() << ") matrix" );
	BOOST_CHECK_MESSAGE( cursor.at_end_of_current_row(), "row(" << cursor.current_row() << ") is not the end of a (" << m.rows() << "," << m.columns() << ") matrix" );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixCursor_start_next_column, MatrixType, matrix_type_list )
{
	using namespace mla::matrix;

	MatrixType m(3,3);
	m.setZero();
	m(0,0) = 1;
	m(1,1) = 2;
	m(2,2) = 3;

	typename MatrixType::Cursor cursor(m);

	cursor.reset();
	cursor.start_next_column();
	BOOST_CHECK_EQUAL( cursor.current_column(), 1 );
	BOOST_CHECK_EQUAL( cursor.current_row(), 0 );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixCursor_start_next_row, MatrixType, matrix_type_list )
{
	using namespace mla::matrix;

	MatrixType m(3,3);
	m.setZero();
	m(0,0) = 1;
	m(1,1) = 2;
	m(2,2) = 3;

	typename MatrixType::Cursor cursor(m);

	cursor.reset();
	cursor.start_next_row();
	BOOST_CHECK_EQUAL( cursor.current_column(), 0 );
	BOOST_CHECK_EQUAL( cursor.current_row(), 1 );
}


BOOST_AUTO_TEST_SUITE_END()

