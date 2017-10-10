#ifndef MLA_SOLVERS_CHOLESKY_HPP
#define MLA_SOLVERS_CHOLESKY_HPP


#include <cmath>

#include <solvers/SolverReturnCodes.h++>
#include <solvers/substitution.h++>

#include <mla/Matrix.h++>
#include <mla/output.h++>

namespace mla
{

/**
Generic cholesky decomposition for dense matrices, works on all matrix types
**/
template<typename Scalar, template<typename> class MatrixStoragePolicy, template<typename> class LMatrixStoragePolicy, template<typename> class VectorStoragePolicy>
ReturnCode 
cholesky(Matrix<Scalar, MatrixStoragePolicy> &A, Vector<Scalar, VectorStoragePolicy> &x, Vector<Scalar, VectorStoragePolicy> &b, Matrix<Scalar, LMatrixStoragePolicy> &L)
{
	if( !A.isSquare() )
	{
		throw LAException("A must be a square matrix");
	}
	if(A.columns() == b.size())
	{
		throw LAException("A.columns() != b.size()");
	}

	Vector<Scalar> S;
	S.resize(A.columns());

	// A=LL^t, factor L
	for(size_t j = 0; j < A.columns(); j++)
	{
		for(size_t i = j; i < A.rows(); i++)
		{
			Scalar Si;
			Si = A.getValue(i,j);
			for(size_t k = 0; k < j; k++)
			{
				Si -= L.getValue(i,k)*L.getValue(j,k);
			}
			S(i) = Si;
		}

		
		Scalar diag = L(j,j) = sqrt(S(j));

		for(size_t i = j+1; i < A.rows(); i++)
		{
			L(i,j) = S(i)/diag;
		}
	}
	// */


	ReturnCode code;
	// Ly = b
	code = forward_substitution(L,x,b);
	if(code != OK)
		return code;

	//Lx=y
	code = back_substitution(L,x,x);
	if(code != OK)
		return code;


	return OK;
}


/**
Generic cholesky decomposition for dense matrices, works on all matrix types
**/
template<typename Scalar, template<typename> class MatrixStoragePolicy, template<typename> class VectorStoragePolicy>
ReturnCode 
cholesky(Matrix<Scalar, MatrixStoragePolicy> &A, Vector<Scalar, VectorStoragePolicy> &x, Vector<Scalar, VectorStoragePolicy> &b, Matrix<Scalar, SparseCRS> &L)
{
	if( !A.isSquare() )
	{
		throw LAException("A must be a square matrix");
	}
	if(A.columns() == b.size())
	{
		throw LAException("A.columns() != b.size()");
	}

	Vector<Scalar, SparseCS> S;
	S.reserve(A.columns());
	S.resize(A.columns());

	size_t *p;
	Scalar *lp;
	size_t *q;
	Scalar *lq;

	// A=LL^t, factor L
	for(size_t j = 0; j < A.columns(); j++)
	{
		S.clear();
		Scalar Si;	// temp variable to avoid calling S(i)
		// unrolling the for(i = j;...) part of the for loop
		Si = A.getValue(j,j);
		p = &L.data.column_index[L.data.row_pointer[j]] ;
		lp = &L.data.values[L.data.row_pointer[j]];
		while(*p < j)
		{
			Si -= (*lp)*(*lp);
	
			p++;
			lp++;
		}
		//S(j) = Si;
		S.push_back(j,Si);

		// continuing the remaining for loop
		for(size_t i = j+1; i < A.rows(); i++)
		{
			Si = A.getValue(i,j);
		
			p = &L.data.column_index[L.data.row_pointer[i]] ;
			lp = &L.data.values[L.data.row_pointer[i]];
			q = &L.data.column_index[L.data.row_pointer[j]] ;
			lq = &L.data.values[L.data.row_pointer[j]] ;

			while( (*p < j) && (*q < j) )
			{
				if(*p == *q)
				{
					Si -= (*lp)*(*lq);
			
					p++;
					lp++;

					q++;
					lq++;
				}
				else if(*p < *q)
				{
					p++;
					lp++;
				}
				else
				{
					q++;
					lq++;
				}
			}
			// */

			//S(i) = Si;
			S.push_back(i,Si);
		}

		Scalar diag;
		diag = L(j,j) = sqrt(S(j));

		for(size_t i = 1; i < S.data.values.size(); i++)
		{
			L(S.data.column_index[i],j) = S.data.values[i]/diag;
		}
	}

	// */


	ReturnCode code;
	// Ly = b
	code = forward_substitution(L,x,b);
	if(code != OK)
		return code;

	//Lx=y
	code = back_substitution(L,x,x);
	if(code != OK)
		return code;

	return OK;
}


/**
Generic cholesky decomposition for dense matrices, works on all matrix types
**/
template<typename Scalar, template<typename> class MatrixStoragePolicy, template<typename> class VectorStoragePolicy>
ReturnCode 
cholesky(Matrix<Scalar, MatrixStoragePolicy> &A, Vector<Scalar, VectorStoragePolicy> &x, Vector<Scalar, VectorStoragePolicy> &b, Matrix<Scalar> &L)
{
	if(A.rows() != A.columns())
	{
		throw LAException("A.rows() != A.columns()");
	}
	if(A.columns() == b.size())
	{
		throw LAException("A.columns() != b.size()");
	}

	Vector<Scalar> S;
	S.resize(A.columns());

	// A=LL^t, factor L
	for(size_t j = 0; j < A.columns(); j++)
	{
		for(size_t i = j; i < A.rows(); i++)
		{
			Scalar Si;
			Scalar *li, *lj;

			Si = A.getValue(i,j);
			li = &L.data.data[i*L.data.t_rows];
			lj = &L.data.data[j*L.data.t_rows];
			for(size_t k = 0; k < j; k++)
			{
				//Si -= L.getValue(i,k)*L.getValue(j,k);
				Si -= (*li)*(*lj);
				++li;
				++lj;
			}
			S(i) = Si;
		}

		
		Scalar diag = L(j,j) = sqrt(S(j));

		for(size_t i = j+1; i < A.rows(); i++)
		{
			L(i,j) = S(i)/diag;
		}
	}
	// */


	ReturnCode code;
	// Ly = b
	code = forward_substitution(L,x,b);
	if(code != OK)
		return code;

	//Lx=y
	code = back_substitution(L,x,x);
	if(code != OK)
		return code;


	return OK;
}


}	// namespace mla

#endif

