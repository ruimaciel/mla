#ifndef MLA_SOLVERS_CONJUGATE_GRADIENT_HPP
#define MLA_SOLVERS_CONJUGATE_GRADIENT_HPP

#include <string>

#include <boost/lexical_cast.hpp>

#include <mla/matrix/all.h++>
#include <mla/vector/all.h++>

#include <mla/operations/level1/dot.h++>
#include <mla/operations/level1/axpy.h++>
#include <mla/operations/level1/scale.h++>
#include <mla/operations/level2/gemv.h++>

#include <mla/solvers/SolverReturnCodes.h++>

#include <mla/output.h++>

namespace mla
{


/**
 * Conjugate gradient method algorithm
 *
 *@param A	Matrix
 *@param x	unknown vector
 *@param b	vector
 *@param delta	error tolerance
 *@param max_iterations	maximum number of iterations allowed
 *
 *@return 
 */
template<typename Scalar, template<typename> class MatrixStoragePolicy, template<typename> class VectorStoragePolicyX, template<typename> class VectorStoragePolicyB >
ReturnCode 
cg(MatrixStoragePolicy<Scalar> &A, VectorStoragePolicyX<Scalar> &x, VectorStoragePolicyB<Scalar> &b, const double delta, int max_iterations) 
{
	if( !A.isSquare() )
	{
		throw LAException("cg: A must be a square matrix");
	}

	if(A.columns() != b.size())
	{
		throw LAException("cg: A.columns() != b.size()");
	}

	x.resize( b.size() );

	VectorStoragePolicyB<Scalar> r, p(A.columns());
	VectorStoragePolicyB<Scalar> Ap(A.columns());
	Scalar dotrr, dotrrnew, alpha;

	//r = b - A*x;
	r = b;
	mla::gemv( (Scalar)(-1.0f), A, x, (Scalar)(1.0f), r );
	p = r;

	dotrr = dot(r,r);

	for (int iter = 0; iter < max_iterations; iter++)
	{
		//Ap = A*p;
		mla::gemv( (Scalar)1.0f, A, p, (Scalar)1.0f, Ap);

		alpha = dotrr/dot(p,Ap);
		//x = x + alpha*p;
		mla::axpy(alpha, p, x);
		
		// r_{k+1] = r_{k} - a*A*p
		mla::axpy( (Scalar)-1.0f*alpha, Ap, r);

		dotrrnew = dot(r,r);

		if(dotrrnew < delta)
		{
			return OK;
		}

		//p = r + (dotrrnew/dotrr)*p;
		mla::scale( dotrrnew/dotrr, p);
		mla::axpy( (Scalar)1.0f, r, p);

		dotrr = dotrrnew;
	}

	return ERR_EXCESSIVE_ITERATIONS;
}


}

#endif
