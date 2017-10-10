#ifndef MLA_MATRIX_STORAGE_POLICY_STATIC_DENSE_ROW_MAJOR_HPP
#define MLA_MATRIX_STORAGE_POLICY_STATIC_DENSE_ROW_MAJOR_HPP

#include <vector>

#include <mla/matrix/traits.h++>
#include <mla/LAException.h++>
#include <mla/vector/SparseCS.h++>



namespace mla
{
namespace matrix
{

/**
StaticDenseRowMajor: a storage policy class intended to implement the interfaces 
specific to the dense, dynamically non-resizeable matrix to be used in the 
Matrix host class.
**/
template<typename Scalar, size_t static_rows, size_t static_columns>
class StaticDenseRowMajor
{
public:
	typedef Scalar scalar_type;

	// type used in row and column vectors
	typedef vector::SparseCS<Scalar>	VectorType;

	std::vector< Scalar >	data;	

public:
	StaticDenseRowMajor();

	/**
	Sets all values to zero
	**/
	void setZero();

	size_t rows() const		{ return static_rows; };
	size_t columns() const		{ return static_columns; };

	/*
	Returns the value in [row,column] 
	*/
	Scalar getValue(size_t row, size_t column) const;

	/**
	 * Sets the value of the element in (row, column)
	 *@param row	the element row
	 *@param column	the element column
	 */
	void setValue(size_t row, size_t column, Scalar value);

	/**
	 * Returns a reference to the element in (row, column)
	 *@param row	the element row
	 *@param column	the element column
	 */
	Scalar & operator() (size_t row, size_t column);


	void resize(size_t row, size_t column);


	/**
	 * Check if matrix is square
	 **/
	bool isSquare() const { return rows() == columns(); }
};



template<typename Scalar, size_t static_rows, size_t static_columns>
StaticDenseRowMajor<Scalar, static_rows, static_columns>::StaticDenseRowMajor()
{
	data.resize(static_rows*static_columns);
}


template<typename Scalar, size_t static_rows, size_t static_columns>
void
StaticDenseRowMajor<Scalar, static_rows, static_columns>::setZero()
{
	using namespace std;

	for(auto i = data.begin(); i != data.end(); i++)
	{
		*i = 0;
	}
}


template<typename Scalar, size_t static_rows, size_t static_columns>
Scalar
StaticDenseRowMajor<Scalar, static_rows, static_columns>::getValue(size_t row, size_t column) const
{
	if(row >= this->rows())
	{
		throw LAException("StaticDenseRowMajor::getValue: row < this->rows()");
	}
	if(column >= this->columns())
	{
		throw LAException("StaticDenseRowMajor::getValue: column < this->columns()");
	}

	return this->data[row*this->rows() + column];
}



template<typename Scalar, size_t static_rows, size_t static_columns>
void
StaticDenseRowMajor<Scalar, static_rows, static_columns>::setValue(size_t row, size_t column, Scalar value)
{
	if(row >= this->rows())
	{
		throw LAException("StaticDenseRowMajor::setValue() row >= this->rows()");
	}
	if(column >= this->columns())
	{
		throw LAException("StaticDenseRowMajor::setValue() column >= this->columns()");
	}

	this->data[row*this->rows() + column] = value;
}



template<typename Scalar, size_t static_rows, size_t static_columns>
Scalar
& StaticDenseRowMajor<Scalar, static_rows, static_columns>::operator() (size_t row, size_t column)
{
	if(row >= this->rows())
	{
		throw LAException("StaticDenseRowMajor::operator() row < this->rows()");
	}
	if(column >= this->columns())
	{
		throw LAException("StaticDenseRowMajor::operator() column < this->columns()");
	}

	return this->data[row*this->rows() + column];
}


template<typename Scalar, size_t static_rows, size_t static_columns>
void
StaticDenseRowMajor<Scalar, static_rows, static_columns>::resize(size_t rows, size_t columns)
{
	static_assert(true, "StaticDenseRowMajor can't be resized.");
}


// Specify traits

template<typename Scalar, size_t static_rows, size_t static_columns>
struct Traits<StaticDenseRowMajor<Scalar, static_rows, static_columns> >
{
	static constexpr bool is_writeable() noexcept { return true; }
	static constexpr bool is_resizeable() noexcept { return false; }
	static constexpr bool is_sparse() noexcept { return false; }
};


}	// namespace matrix
}	// namespace mla
#endif
