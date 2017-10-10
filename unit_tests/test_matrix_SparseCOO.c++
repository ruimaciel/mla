#define BOOST_TEST_MODULE matrix

#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


#include <mla/matrix/SparseCOO.h++>


typedef boost::mpl::list<
	float,
	double
> scalar_list;


template <typename Scalar>
using MatrixType = mla::matrix::SparseCOO<Scalar>;


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


BOOST_AUTO_TEST_CASE_TEMPLATE( splitOutBottom_diagonal_matrix, Scalar, scalar_list )
{
	MatrixType<Scalar> m(5,5);
	m.setZero();

	m.setValue(0,0, 0.0);	
	m.setValue(1,1, 1.0);	
	m.setValue(2,2, 2.0);	
	m.setValue(3,3, 3.0);	
	m.setValue(4,4, 4.0);	

	MatrixType<Scalar> m2 = m.splitOutBottom(2);

	// test NNZ size
	BOOST_CHECK_EQUAL(m.nnz(), 2);
	BOOST_CHECK_EQUAL(m2.nnz(), 3);

	// test matrix size
	BOOST_CHECK_EQUAL(m.rows(), 2);
	BOOST_CHECK_EQUAL(m.columns(), 5);
	BOOST_CHECK_EQUAL(m2.rows(), 3);
	BOOST_CHECK_EQUAL(m2.columns(), 5);

	// test matrix elements
	BOOST_CHECK_EQUAL(m.getValue(0, 0), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(1, 1), 1.0);
	BOOST_CHECK_EQUAL(m2.getValue(0, 2), 2.0);
	BOOST_CHECK_EQUAL(m2.getValue(1, 3), 3.0);
	BOOST_CHECK_EQUAL(m2.getValue(2, 4), 4.0);
}


BOOST_AUTO_TEST_CASE_TEMPLATE( splitOutRight_diagonal_matrix, Scalar, scalar_list )
{
	MatrixType<Scalar> m(5,5);
	m.setZero();

	m.setValue(0,0, 0.0);	
	m.setValue(1,1, 1.0);	
	m.setValue(2,2, 2.0);	
	m.setValue(3,3, 3.0);	
	m.setValue(4,4, 4.0);	

	MatrixType<Scalar> m2 = m.splitOutRight(2);

	// test NNZ size
	BOOST_CHECK_EQUAL(m.nnz(), 2);
	BOOST_CHECK_EQUAL(m2.nnz(), 3);

	// test matrix size
	BOOST_CHECK_EQUAL(m.rows(), 5);
	BOOST_CHECK_EQUAL(m.columns(), 2);
	BOOST_CHECK_EQUAL(m2.rows(), 5);
	BOOST_CHECK_EQUAL(m2.columns(), 3);

	// test matrix elements
	BOOST_CHECK_EQUAL(m.getValue(0, 0), 0.0);
	BOOST_CHECK_EQUAL(m.getValue(1, 1), 1.0);
	BOOST_CHECK_EQUAL(m2.getValue(2, 0), 2.0);
	BOOST_CHECK_EQUAL(m2.getValue(3, 1), 3.0);
	BOOST_CHECK_EQUAL(m2.getValue(4, 2), 4.0);
}




BOOST_AUTO_TEST_SUITE_END()

