#define BOOST_TEST_MODULE matrix

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/all.h++>
#include <mla/matrix/Diagonal.h++>


typedef boost::mpl::list<
	float,
	double
> scalar_list;


template <typename Scalar>
using MatrixType = mla::matrix::Diagonal<Scalar>;


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


BOOST_AUTO_TEST_SUITE_END()

