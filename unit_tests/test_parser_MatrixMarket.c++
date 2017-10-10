#define BOOST_TEST_MODULE parsers

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <mla/matrix/all.h++>
#include <mla/parsers/MatrixMarket.h++>

#include "config_paths.h"

typedef boost::mpl::list<
	mla::matrix::DenseRowMajor<float>,
	mla::matrix::DenseRowMajor<double>,
	mla::matrix::SparseDOK<float>,
	mla::matrix::SparseDOK<double>
> matrix_type_list;


typedef boost::mpl::list<
//	float,
	double
> scalar_list;


template <typename Scalar>
using MatrixType = mla::matrix::SparseDOK<Scalar>;


BOOST_AUTO_TEST_SUITE(parser_MatrixMarket)


BOOST_AUTO_TEST_CASE_TEMPLATE( test_parser_coordinate_real_general, Scalar, scalar_list )
{
	using namespace mla;

	MatrixMarket parser;
	std::stringstream ss;
	ss.str("%%MatrixMarket matrix coordinate real general\n\
5 5 2\n\
1 2 3.0\n\
4 5 6.0\n\
");


	std::vector<MatrixMarket::Symbol> expected_token_vector = { 
		MatrixMarket::TS_COMMENT,
		MatrixMarket::TS_COMMENT,
		MatrixMarket::TS_STRING_MATRIX_MARKET,
		MatrixMarket::TS_STRING_MATRIX,
		MatrixMarket::TS_STRING_COORDINATE,
		MatrixMarket::TS_STRING_REAL,
		MatrixMarket::TS_STRING_GENERAL,
		MatrixMarket::TS_EOL
	};

	// test the tokens
	for(MatrixMarket::Symbol expected_token: expected_token_vector)
	{
		MatrixMarket::Symbol received_token = parser.lexer_header(ss);
		BOOST_CHECK_EQUAL( received_token, expected_token);
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE( test_lexer_body, Scalar, scalar_list )
{
	using namespace mla;

	MatrixMarket parser;
	std::stringstream ss;
	ss.str("0\n1\n%some crap\n1.0\n-1\n1.234e-56\n");


	std::vector<MatrixMarket::Symbol> expected_token_vector = { 
		MatrixMarket::TS_INTEGER_NUMBER,
		MatrixMarket::TS_EOL,
		MatrixMarket::TS_INTEGER_NUMBER,
		MatrixMarket::TS_EOL,
		MatrixMarket::TS_COMMENT,
		MatrixMarket::TS_EOL,
		MatrixMarket::TS_REAL_NUMBER,
		MatrixMarket::TS_EOL,
		MatrixMarket::TS_REAL_NUMBER,
		MatrixMarket::TS_EOL,
		MatrixMarket::TS_REAL_NUMBER,
		MatrixMarket::TS_EOL
	};

	// test the tokens
	for(MatrixMarket::Symbol expected_token: expected_token_vector)
	{
		MatrixMarket::Symbol received_token = parser.lexer_body(ss);
		BOOST_CHECK_EQUAL( received_token, expected_token);
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE( test_lexer_body_integer_number, Scalar, scalar_list )
{
	using namespace mla;

	MatrixMarket parser;
	std::stringstream ss;
	ss.str("0 1 2");


	std::vector<unsigned long> expected_integer_number_vector = { 
		0,
		1,
		2,
	};

	// test the tokens
	for(unsigned int expected_number: expected_integer_number_vector)
	{
		parser.lexer_body(ss);
		BOOST_CHECK_EQUAL( expected_number, parser.integer_number);
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE( test_lexer_body_real_number, Scalar, scalar_list )
{
	using namespace mla;

	MatrixMarket parser;
	std::stringstream ss;
	ss.str("0.0 1.0 2.0    3.1   ");


	std::vector<double> expected_integer_number_vector = { 
		0.0,
		1.0,
		2.0,
		3.1
	};

	// test the tokens
	for(double expected_number: expected_integer_number_vector)
	{
		parser.lexer_body(ss);
		BOOST_CHECK_CLOSE( expected_number, parser.real_number, 1e-5);
	}
}


BOOST_AUTO_TEST_CASE_TEMPLATE( basic_coordinate_real_general_matrix, MatrixType, matrix_type_list )
{
	using namespace mla;

	MatrixMarket parser;
	MatrixType matrix;
	std::stringstream ss;

	ss.str("%%MatrixMarket matrix coordinate real general \n\
8 9 2 \n\
1 3 4.0 \n\
5 6 7.0 \n\
");

	parser.parse(ss, matrix);

	// test
	BOOST_CHECK_EQUAL( matrix.rows(), 8 );
	BOOST_CHECK_EQUAL( matrix.columns(), 9 );
	BOOST_CHECK_CLOSE( matrix.getValue(0, 2), 4.0, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(4, 5), 7.0, 1e-5);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( basic_coordinate_real_symmetric_matrix, MatrixType, matrix_type_list )
{
	using namespace mla;

	MatrixMarket parser;
	MatrixType matrix;
	std::stringstream ss;

	ss.str("%%MatrixMarket matrix coordinate real symmetric \n\
9 9 2 \n\
1 3 4.0 \n\
5 6 7.0 \n\
");

	parser.parse(ss, matrix);

	// test
	BOOST_CHECK_EQUAL( matrix.rows(), 9 );
	BOOST_CHECK_EQUAL( matrix.columns(), 9 );
	BOOST_CHECK_CLOSE( matrix.getValue(0, 2), 4.0, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(2, 0), 4.0, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(4, 5), 7.0, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(5, 4), 7.0, 1e-5);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( basic_coordinate_real_skew_symmetric_matrix, MatrixType, matrix_type_list )
{
	using namespace mla;

	MatrixMarket parser;
	MatrixType matrix;
	std::stringstream ss;

	ss.str("%%MatrixMarket matrix coordinate real skew-symmetric \n\
9 9 2 \n\
1 3 4.0 \n\
5 6 7.0 \n\
");

	parser.parse(ss, matrix);

	// test
	BOOST_CHECK_EQUAL( matrix.rows(), 9 );
	BOOST_CHECK_EQUAL( matrix.columns(), 9 );
	BOOST_CHECK_CLOSE( matrix.getValue(0, 2),  4.0, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(2, 0), -4.0, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(4, 5),  7.0, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(5, 4), -7.0, 1e-5);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( parse_coordinate_test01, MatrixType, matrix_type_list )
{
	using namespace mla;

	MatrixMarket parser;
	MatrixType matrix;

	std::string file_path = FULL_MATRIX_MARKET_FILES_PATH "coordinate/test01.mtx";
	
	std::ifstream file(file_path, std::ifstream::in);

	BOOST_CHECK( file.is_open() );

	// parse
	parser.parse(file, matrix);

	BOOST_CHECK_EQUAL( matrix.rows(), 5 );
	BOOST_CHECK_EQUAL( matrix.columns(), 5 );

	BOOST_CHECK_CLOSE( matrix.getValue(0, 0),  1.0, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(1, 1),  10.5, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(3, 1),  250.5, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(2, 2),  0.015, 1e-5);

	BOOST_CHECK_CLOSE( matrix.getValue(0, 3),  6.0, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(3, 3),  -280.0, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(3, 4),  33.32, 1e-5);
	BOOST_CHECK_CLOSE( matrix.getValue(4, 4),  12.0, 1e-5);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( parse_coordinate_bcsstk14, MatrixType, matrix_type_list )
{
	using namespace mla;

	/**
	 * Parse file taken from:
	 * http://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc2/bcsstk14.html
	 **/

	MatrixMarket parser;
	MatrixType matrix;

	std::string file_path = FULL_MATRIX_MARKET_FILES_PATH "coordinate/bcsstk14.mtx";
	
	std::ifstream file(file_path, std::ifstream::in);

	BOOST_CHECK( file.is_open() );

	// parse
	parser.parse(file, matrix);

	BOOST_CHECK_EQUAL( matrix.rows(), 1806 );
	BOOST_CHECK_EQUAL( matrix.columns(), 1806 );
}


BOOST_AUTO_TEST_SUITE_END()

