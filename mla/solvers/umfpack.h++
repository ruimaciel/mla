#ifndef MLA_SOLVERS_UMFPACK_HPP
#define MLA_SOLVERS_UMFPACK_HPP


#include <vector>
#include <suitesparse/umfpack.h>

#include <mla/solvers/SolverReturnCodes.h++>
#include <mla/ProgressIndicatorStrategy.h++>

#include <mla/matrix/SparseCCS.h++>
#include <mla/vector/Dense.h++>


namespace mla
{

/**
Umfpack routine
**/
ReturnCode 
umfpack(matrix::SparseCCS<double> &A, vector::Dense<double> &x, vector::Dense<double> &b, ProgressIndicatorStrategy *)
{
	if( !A.isSquare() )
	{
		throw LAException("A must be a square matrix");
	}
	if(A.columns() != b.size())
	{
		throw LAException("A.columns() != b.size()");
	}


	x.resize(A.rows());


	void *Numeric = NULL;
	void *Symbolic = NULL;

	double *null = (double *) NULL ;

	long int row = A.rows();
	long int column = A.columns();

	std::vector<long int> Ap( A.data.column_pointer.begin(), A.data.column_pointer.end());
	std::vector<long int> Ai( A.data.row_index.begin(), A.data.row_index.end());

	(void) umfpack_dl_symbolic (row, column, &Ap[0], &Ai[0], &A.data.values[0], &Symbolic, null, null) ;
	(void) umfpack_dl_numeric (&Ap[0], &Ai[0], &A.data.values[0], Symbolic, &Numeric, null, null);
	
	umfpack_dl_free_symbolic (&Symbolic) ;

	(void) umfpack_dl_solve (UMFPACK_A, &Ap[0], &Ai[0], &A.data.values[0], &x.data[0], &b.data[0], Numeric, null, null) ;
	umfpack_dl_free_numeric (&Numeric) ;

	#undef Ap
	#undef Ai


	return OK;
}


}	// namespace mla

#endif
