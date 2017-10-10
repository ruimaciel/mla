#define BOOST_TEST_MODULE matrix

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/SparseCRS.h++>


typedef boost::mpl::list<
	float,
	double
> scalar_list;


template <typename Scalar>
using MatrixType = mla::matrix::SparseCRS<Scalar>;


BOOST_AUTO_TEST_SUITE(test_matrix)


BOOST_AUTO_TEST_CASE_TEMPLATE( Declaration, Scalar, scalar_list )
{
	MatrixType<Scalar> m(3,3);

	BOOST_CHECK_EQUAL(m.rows(), 3);
	BOOST_CHECK_EQUAL(m.columns(), 3);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( Clearing, Scalar, scalar_list )
{
	MatrixType<Scalar> m(3,3);
	m.setZero();

	for(unsigned int i = 0; i < 3; i++)
	{
		for(unsigned int j = 0; j < 3; j++)
		{
			BOOST_CHECK_EQUAL(m.getValue(i, j), 0.0);
		}
	}
	
}


BOOST_AUTO_TEST_CASE_TEMPLATE( Assigning, Scalar, scalar_list )
{
	MatrixType<Scalar> m(3,3);
	m.setZero();

	m(1,1) = 1.0;	

	BOOST_CHECK_EQUAL(m.getValue(1, 1), 1.0);

	BOOST_CHECK_EQUAL(m.getValue(0, 1), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(2, 1), 0.0);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( get_row, Scalar, scalar_list )
{
	MatrixType<Scalar> m(3,3);
	m.setZero();

	m(0,0) = 1.0;	
	m(1,1) = 2.0;	
	m(2,2) = 3.0;	


	for(size_t i = 0; i < m.rows(); i++)
	{
		typename MatrixType<Scalar>::VectorType v = m.getRow(i);
		for(size_t j = 0; j < m.columns(); j++)
		{
			BOOST_CHECK_EQUAL(v.getValue(j), m.getValue(i,j));
		}
	}

}


BOOST_AUTO_TEST_CASE_TEMPLATE( Incrementing_21, Scalar, scalar_list )
{
	MatrixType<Scalar> m(3,3);
	m.setZero();

	m(2,1) += 1.0;	

	BOOST_CHECK_EQUAL(m.getValue(0, 0), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(0, 1), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(0, 2), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(1, 0), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(1, 1), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(1, 2), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(2, 0), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(2, 1), 1.0);
	BOOST_CHECK_EQUAL(m.getValue(2, 2), 0.0);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( Incrementing_02, Scalar, scalar_list )
{
	MatrixType<Scalar> m(3,3);
	m.setZero();

	m(0,2) += 1.0;	

	BOOST_CHECK_EQUAL(m.getValue(0, 0), (Scalar)0);
	BOOST_CHECK_EQUAL(m.getValue(0, 1), (Scalar)0);
	BOOST_CHECK_EQUAL(m.getValue(0, 2), (Scalar)1);
	BOOST_CHECK_EQUAL(m.getValue(1, 0), (Scalar)0);
	BOOST_CHECK_EQUAL(m.getValue(1, 1), (Scalar)0);
	BOOST_CHECK_EQUAL(m.getValue(1, 2), (Scalar)0);
	BOOST_CHECK_EQUAL(m.getValue(2, 0), (Scalar)0);
	BOOST_CHECK_EQUAL(m.getValue(2, 1), (Scalar)0);
	BOOST_CHECK_EQUAL(m.getValue(2, 2), (Scalar)0);
}


BOOST_AUTO_TEST_SUITE_END()

