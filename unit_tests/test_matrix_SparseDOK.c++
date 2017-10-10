#define BOOST_TEST_MODULE matrix

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/SparseDOK.h++>


typedef boost::mpl::list<
	float,
	double
> scalar_list;


template <typename Scalar>
using MatrixType = mla::matrix::SparseDOK<Scalar>;


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


BOOST_AUTO_TEST_CASE_TEMPLATE( Permute_rows, Scalar, scalar_list )
{
	MatrixType<Scalar> m(5,5);
	m.setZero();

	m(1,1) = 1.0;	
	m(2,2) = 2.0;	
	m(3,3) = 3.0;	
	m(4,4) = 4.0;	

	m.permuteRows(2,3);

	BOOST_CHECK_EQUAL(m.getValue(0, 0), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(1, 1), 1.0);
	BOOST_CHECK_EQUAL(m.getValue(2, 2), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(2, 3), 3.0);
	BOOST_CHECK_EQUAL(m.getValue(3, 2), 2.0);
	BOOST_CHECK_EQUAL(m.getValue(3, 3), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(4, 4), 4.0);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( Permute_columns, Scalar, scalar_list )
{
	MatrixType<Scalar> m(5,5);
	m.setZero();

	m(1,1) = 1.0;	
	m(2,2) = 2.0;	
	m(3,3) = 3.0;	
	m(4,4) = 4.0;	

	m.permuteColumns(2,3);

	BOOST_CHECK_EQUAL(m.getValue(0, 0), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(1, 1), 1.0);
	BOOST_CHECK_EQUAL(m.getValue(2, 2), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(2, 3), 2.0);
	BOOST_CHECK_EQUAL(m.getValue(3, 2), 3.0);
	BOOST_CHECK_EQUAL(m.getValue(3, 3), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(4, 4), 4.0);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( symmetric_reorder_identity_matrix, Scalar, scalar_list )
{
	MatrixType<Scalar> m(5,5);
	m.setZero();

	m(1,1) = 1.0;	
	m(2,2) = 2.0;	
	m(3,3) = 3.0;	
	m(4,4) = 4.0;	

	std::vector<size_t> new_order{0,1,2,3,4};

	m.symmetricReorder(new_order);

	BOOST_CHECK_EQUAL(m.getValue(0, 0), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(1, 1), 1.0);
	BOOST_CHECK_EQUAL(m.getValue(2, 2), 2.0);
	BOOST_CHECK_EQUAL(m.getValue(2, 3), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(3, 2), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(3, 3), 3.0);
	BOOST_CHECK_EQUAL(m.getValue(4, 4), 4.0);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( symmetric_reorder_identity_matrix_2, Scalar, scalar_list )
{
	MatrixType<Scalar> m(5,5);
	m.setZero();

	m(1,1) = 1.0;	
	m(2,2) = 2.0;	
	m(3,3) = 3.0;	
	m(4,4) = 4.0;	

	std::vector<size_t> new_order{0,2,1,3,4};

	m.symmetricReorder(new_order);

	BOOST_CHECK_EQUAL(m.getValue(0, 0), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(1, 1), 2.0);
	BOOST_CHECK_EQUAL(m.getValue(2, 2), 1.0);
	BOOST_CHECK_EQUAL(m.getValue(2, 3), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(3, 2), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(3, 3), 3.0);
	BOOST_CHECK_EQUAL(m.getValue(4, 4), 4.0);
}



BOOST_AUTO_TEST_SUITE_END()

