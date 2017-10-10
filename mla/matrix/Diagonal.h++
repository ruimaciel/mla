#ifndef MLA_MATRIX_STORAGE_POLICY_DIAGONAL_HPP
#define MLA_MATRIX_STORAGE_POLICY_DIAGONAL_HPP

#include <vector>

#include <mla/matrix/traits.h++>
#include <mla/LAException.h++>
#include <mla/MatrixCursor.h++>
#include <mla/vector/SparseCS.h++>


namespace mla
{
namespace matrix
{

/**
Diagonal: a storage policy class for the Matrix host class.
This class implements the interface for the diagonal matrix format
**/
template<typename Scalar>
class Diagonal
{
public:
	typedef Scalar scalar_type;

	// type used in row and column vectors
	typedef vector::SparseCS<Scalar>	VectorType;

	struct Data
	{
		size_t	t_rows;		// number of rows
		size_t	t_columns;	// number of columns

		std::vector<Scalar> data;
	} data;

public:
	Diagonal(size_t rows = 0, size_t columns = 0);

	/**
	Sets all values to zero
	**/
	void setZero();

	size_t rows() const		{ return data.t_rows; };
	size_t columns() const		{ return data.t_columns; };

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
	Sets ones on the main diagonal and zeros elsewhere
	**/
	void setEye();


	/**
	 * Check if matrix is square
	 **/
	bool isSquare() const { return rows() == columns(); }


	/**
	 * An element iterator class that is a convenient generic method
	 * to iterate through matrix elements
	 */
	class Cursor
		: public MatrixCursor<Scalar, Diagonal>
	{
	public:
		Cursor(Diagonal &m);

		/**
		 * resets the cursor to the first non-null term of the matrix
		 */
		void reset() override final;

		/**
		 * returns the element
		 **/
		Scalar element() override final;

		/**
		 * Returns the column currently being pointed out by the cursor
		 */
		size_t current_column() const override final;

		/**
		 * Returns the row currently being pointed out by the cursor
		 */
		size_t current_row() const override final;

		/**
		 * Check if cursor is pointing to the first element of the current column
		 */
		bool at_beginning_of_column() const override final;
		bool at_beginning_of_row() const override final;

		/**
		 * Check if cursor is pointing to the last element of the current column
		 */
		bool at_end_of_current_column() const override final;
		bool at_end_of_current_row() const override final;

		/**
		 * Check if cursor is pointing to an element pointing to the last column
		 */
		bool at_end_of_columns() const override final;
		bool at_end_of_rows() const override final;


		/**
		 * resets the cursor to the first non-null term of the matrix
		 */
		void start_current_column() override final;
		void start_current_row() override final;

		/**
		 * resets the cursor to the first non-null term of the matrix
		 */
		void start_current_column_nn() override final;
		void start_current_row_nn() override final;

		/**
		 * Moves cursor to the first non-null element of the subsequent column
		 */
		void start_next_column_nn() override final;
		void start_next_row_nn() override final;

		/**
		 * Moves cursor to next non-null element in the same row
		 */
		void increment_column() override final;
		void increment_row() override final;

	};


	/**
	 * returns a MatrixCursor
	 */
	Cursor cursor()	{	return	typename Diagonal<Scalar>::Cursor(*this);	}

};



template<typename Scalar>
Diagonal<Scalar>::Diagonal(size_t rows, size_t columns)
{
	resize(rows, columns);
}



template<typename Scalar>
void
Diagonal<Scalar>::setZero()
{
	using namespace std;
	for(typename std::vector<Scalar>::iterator i = data.data.begin(); i != data.data.end(); i++)
	{
		*i = 0;
	}
}


template<typename Scalar>
Scalar
Diagonal<Scalar>::getValue(size_t row, size_t column) const
{
	using namespace std;

	if(row >= data.t_rows)
	{
		throw LAException("row < data.t_rows");
	}
	if(column >= data.t_columns)
	{
		throw LAException("column < data.t_columns");
	}


	if(row != column)
		return 0;

	return data.data[row];
}


template<typename Scalar>
void
Diagonal<Scalar>::setValue(size_t row, size_t column, Scalar value)
{
	if(row >= data.t_rows)
	{
		throw LAException("Diagonal::setValue(): row >= data.t_rows");
	}
	if(column >= data.t_columns)
	{
		throw LAException("Diagonal::setValue(): columns >= data.t_columns");
	}
	if(row != column)
	{
		throw LAException("Diagonal: cannot set values in non-diagonal elements of a diagonal matrix");
	}

	this->data.data[row] = value;
}


template<typename Scalar>
Scalar &
Diagonal<Scalar>::operator() (size_t row, size_t column)
{
	if(row >= data.t_rows)
	{
		throw LAException("row < data.t_rows");
	}
	if(column >= data.t_columns)
	{
		throw LAException("column < data.t_columns");
	}
	if(row != column)
	{
		throw LAException("cannot set values in non-diagonal elements of a diagonal matrix");
	}

	return this->data.data[row];
}


template<typename Scalar>
void
Diagonal<Scalar>::resize(size_t rows, size_t columns)
{
	using namespace std;


	rows < columns ? data.data.resize(rows): data.data.resize(columns);

	data.t_rows = rows;
	data.t_columns = columns;
}


template<typename Scalar>
void
Diagonal<Scalar>::setEye()
{
	for(auto &e: data.data)
	{
		e = (Scalar)1;
	}
}


/**
 * Start of cursor definition
 **/

template<typename Scalar>
Diagonal<Scalar>::Cursor::Cursor(Diagonal &m)
	: MatrixCursor<Scalar, Diagonal> (m)
{
	this->reset();
}


template<typename Scalar>
void
Diagonal<Scalar>::Cursor::start_current_column()
{
	this->m_j = 0;
}


template<typename Scalar>
void
Diagonal<Scalar>::Cursor::start_current_row()
{
	this->m_i = 0;
}


template<typename Scalar>
void
Diagonal<Scalar>::Cursor::start_current_column_nn()
{
	this->m_i = this->m_j;
}


template<typename Scalar>
void
Diagonal<Scalar>::Cursor::start_current_row_nn()
{
	this->m_j = this->m_i;
}


template<typename Scalar>
void
Diagonal<Scalar>::Cursor::reset()
{
	this->m_i = 0;
	this->m_j = 0;
}


template<typename Scalar>
Scalar
Diagonal<Scalar>::Cursor::element()
{
	if( this->m_i != this->m_j )
	{
		return 0.0;
	}
	else
	{
		return this->m_matrix(this->m_i, this->m_j);
	}
}


template<typename Scalar>
size_t
Diagonal<Scalar>::Cursor::current_row() const
{
	return this->m_i;
}


template<typename Scalar>
size_t
Diagonal<Scalar>::Cursor::current_column() const
{
	return this->m_j;
}



template<typename Scalar>
bool
Diagonal<Scalar>::Cursor::at_beginning_of_column() const
{
	return this->m_j == 0;
}

	
template<typename Scalar>
bool
Diagonal<Scalar>::Cursor::at_beginning_of_row() const
{
	return this->m_j == 0;
}

template<typename Scalar>
bool
Diagonal<Scalar>::Cursor::at_end_of_current_column() const
{
	return !(this->m_j < this->m_matrix.columns() );
}

	
template<typename Scalar>
bool
Diagonal<Scalar>::Cursor::at_end_of_current_row() const
{
	return !(this->m_i < this->m_matrix.rows() );
}


template<typename Scalar>
bool
Diagonal<Scalar>::Cursor::at_end_of_columns() const
{
	return !(this->m_j < this->m_matrix.columns() );
}

	
template<typename Scalar>
bool
Diagonal<Scalar>::Cursor::at_end_of_rows() const
{
	return !(this->m_i < this->m_matrix.rows() );
}


template<typename Scalar>
void
Diagonal<Scalar>::Cursor::start_next_column_nn()
{
	this->m_j++;
	this->m_i = this->m_j;
}


template<typename Scalar>
void
Diagonal<Scalar>::Cursor::start_next_row_nn()
{
	this->m_i++;
	this->m_j = this->m_i;
}


template<typename Scalar>
void
Diagonal<Scalar>::Cursor::increment_column()
{
	this->m_j++; 
}


template<typename Scalar>
void
Diagonal<Scalar>::Cursor::increment_row()
{
	this->m_i++; 
}


// Specify traits

template<typename Scalar>
struct Traits<Diagonal<Scalar> >
{
	static constexpr bool is_writeable() noexcept { return false; }
	static constexpr bool is_resizeable() noexcept { return true; }
	static constexpr bool is_sparse() noexcept { return false; }
};


}	// namespace matrix
}	// namespace mla

#endif
