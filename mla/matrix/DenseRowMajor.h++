#ifndef MLA_MATRIX_STORAGE_POLICY_DENSE_ROW_MAJOR_HPP
#define MLA_MATRIX_STORAGE_POLICY_DENSE_ROW_MAJOR_HPP

#include <vector>

#include <mla/matrix/traits.h++>
#include <mla/LAException.h++>
#include <mla/MatrixCursor.h++>
#include <mla/matrix/RowCursor.h++>
#include <mla/vector/SparseCS.h++>


namespace mla
{
namespace matrix
{

/**
DenseRowMajor: a storage policy class intended to implement the interfaces 
specific to the dense, row-major ordering, dynamically resizeable matrix to 
be used in the Matrix host class.
**/
template<typename Scalar>
class DenseRowMajor
{
public:
	typedef Scalar scalar_type;

	// type used in row and column vectors
	typedef vector::SparseCS<Scalar>	VectorType;

	struct Data
	{
		size_t	n_rows;		// number of rows
		size_t	n_columns;	// number of columns

		std::vector< Scalar >	element_vector;	

		size_t getElementIndex(size_t i, size_t j) const {return i+ n_rows*j;}
	} data;

public:
	DenseRowMajor(size_t rows = 0, size_t columns = 0);

	/**
	Sets all values to zero
	**/
	void setZero();

	size_t rows() const		{ return data.n_rows; };
	size_t columns() const		{ return data.n_columns; };

	/*
	Returns the value in [row,column] 
	*/
	Scalar getValue(size_t i, size_t j) const;


	/**
	 * Returns the row vector
	 **/
	VectorType getRow(size_t i) const;

	/**
	 * Returns the column vector
	 **/
	VectorType getColumn(size_t j) const;

	/**
	 * Sets the value of the element in (i,j)
	 *@param i	the element row
	 *@param j	the element column
	 */
	void setValue(size_t i, size_t j, Scalar value);

	/**
	Sets ones on the main diagonal and zeros elsewhere
	**/
	void setEye();

	/**
	 * Returns a reference to the element in (row, column)
	 *@param i	the element row
	 *@param j	the element column
	 */
	Scalar & operator() (size_t i, size_t j);


	void resize(size_t row, size_t column);


	/**
	 * Check if matrix is square
	 **/
	bool isSquare() const { return rows() == columns(); }


	/**
	 * An element iterator class that is a convenient generic method
	 * to iterate through matrix elements
	 */
	class Cursor
		: public MatrixCursor<Scalar, DenseRowMajor>
	{
	public:
		Cursor(DenseRowMajor &m);

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
	Cursor	cursor()	{	return	typename DenseRowMajor<Scalar>::Cursor(*this);	}


	/**
	 * A cursor design pattern for row-wise iterations
	 **/
	class RowCursor 
		: public matrix::RowCursor<Scalar>
	{
	protected:
		DenseRowMajor &m_matrix;
		size_t m_i, m_j;

	public:
		RowCursor(DenseRowMajor &m);

		void reset() override;
		Scalar getElement() const override;
		void setElement(Scalar value) const override;

		size_t getCurrentRow() const override	{ return this->m_i; };

		/**
		 * Returns the column currently being pointed out by the cursor
		 */
		size_t getCurrentColumn() const override { return this->m_j; };

		/**
		 * Moves cursor to next non-null element in the same row
		 */
		void getNextNonNullRowElement() override;
		bool hasNextNonNullRowElement() const override;

		void startNextNonNullRow() override;
		bool hasNextNonNullRow() const override;
	};

	/**
	 * returns a MatrixCursor
	 */
	RowCursor	getRowCursor()	{	return	typename DenseRowMajor<Scalar>::RowCursor(*this);	}
};



template<typename Scalar>
DenseRowMajor<Scalar>::DenseRowMajor(size_t rows, size_t columns)
{
	resize(rows,columns);
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::setZero()
{
	using namespace std;

	data.element_vector.clear();
	size_t n_elements = rows()*columns();
	data.element_vector.resize(n_elements);
}


template<typename Scalar>
Scalar
DenseRowMajor<Scalar>::getValue(size_t i, size_t j) const
{
	if(i >= this->rows())
	{
		throw LAException("DenseRowMajor::getValue: i < rows()");
	}
	if(j >= this->columns())
	{
		throw LAException("DenseRowMajor::getValue: j < columns()");
	}

	return this->data.element_vector[ data.getElementIndex(i,j) ];
}


template<typename Scalar>
typename DenseRowMajor<Scalar>::VectorType 
DenseRowMajor<Scalar>::getRow(size_t i) const
{
	if(i >= this->rows() )
	{
		throw LAException("DenseRowMajor:: i > rows()");
	}


	DenseRowMajor<Scalar>::VectorType row_vector(this->columns());

	for(size_t j = 0; j < this->columns(); j++)
	{
		row_vector.setValue(j, this->getValue(i, j) );
	}

	return row_vector;
}


template<typename Scalar>
typename DenseRowMajor<Scalar>::VectorType 
DenseRowMajor<Scalar>::getColumn(size_t j) const
{
	if(j >= this->columns() )
	{
		throw LAException("DenseRowMajor:: j > this.columns()");
	}


	DenseRowMajor<Scalar>::VectorType column_vector(this->rows());

	for(size_t i = 0; i < this->rows(); i++)
	{
		column_vector.setValue(i, this->getValue(i, j) );
	}

	return column_vector;
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::setValue(size_t i, size_t j, Scalar value)
{
	if(i >= this->rows())
	{
		throw LAException("DenseRowMajor::setValue() i < rows()");
	}
	if(j >= this->columns())
	{
		throw LAException("DenseRowMajor::setValue() j < columns()");
	}

	this->data.element_vector[data.getElementIndex(i,j)] = value;
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::setEye()
{
	for(size_t i = 0; i < this->rows(); i++)
	{
		for(size_t j = 0; j < this->columns(); j++)
		{
			this->setValue(i,j, i==j? (Scalar)1: (Scalar)0);
		}
	}
}



template<typename Scalar>
Scalar & 
DenseRowMajor<Scalar>::operator() (size_t i, size_t j)
{
	if(i >= this->rows())
	{
		throw LAException("DenseRowMajor::operator() i >= rows()");
	}
	if(j >= this->columns())
	{
		throw LAException("DenseRowMajor::operator() j >= columns()");
	}

	return this->data.element_vector[data.getElementIndex(i,j)];
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::resize(size_t rows, size_t columns)
{
	using namespace std;

	data.n_rows = rows;
	data.n_columns = columns;

	size_t n_elements = this->rows()*this->columns();
	data.element_vector.resize(n_elements);
}



template<typename Scalar>
DenseRowMajor<Scalar>::Cursor::Cursor(DenseRowMajor &m)
	: MatrixCursor<Scalar, DenseRowMajor> (m)
{
	this->reset();
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::Cursor::start_current_column()
{
	this->m_j = 0;
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::Cursor::start_current_row()
{
	this->m_i = 0;
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::Cursor::start_current_column_nn()
{
	this->start_current_column();
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::Cursor::start_current_row_nn()
{
	this->start_current_row();
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::Cursor::reset()
{
	this->m_i = 0;
	this->m_j = 0;
}


template<typename Scalar>
Scalar
DenseRowMajor<Scalar>::Cursor::element()
{
	return this->m_matrix(this->m_i, this->m_j);
}


template<typename Scalar>
size_t
DenseRowMajor<Scalar>::Cursor::current_row() const
{
	return this->m_i;
}


template<typename Scalar>
size_t
DenseRowMajor<Scalar>::Cursor::current_column() const
{
	return this->m_j;
}



template<typename Scalar>
bool
DenseRowMajor<Scalar>::Cursor::at_beginning_of_column() const
{
	return this->m_j == 0;
}

	
template<typename Scalar>
bool
DenseRowMajor<Scalar>::Cursor::at_beginning_of_row() const
{
	return this->m_j == 0;
}

template<typename Scalar>
bool
DenseRowMajor<Scalar>::Cursor::at_end_of_current_column() const
{
	return !(this->m_i < this->m_matrix.rows() );
}

	
template<typename Scalar>
bool
DenseRowMajor<Scalar>::Cursor::at_end_of_current_row() const
{
	return !(this->m_j < this->m_matrix.columns() );
}


template<typename Scalar>
bool
DenseRowMajor<Scalar>::Cursor::at_end_of_columns() const
{
	return !(this->m_j < this->m_matrix.columns() );
}

	
template<typename Scalar>
bool
DenseRowMajor<Scalar>::Cursor::at_end_of_rows() const
{
	return !(this->m_i < this->m_matrix.rows() );
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::Cursor::start_next_column_nn()
{
	this->start_next_column();
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::Cursor::start_next_row_nn()
{
	this->start_next_row();
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::Cursor::increment_column()
{
	this->m_j++; 
}


template<typename Scalar>
void
DenseRowMajor<Scalar>::Cursor::increment_row()
{
	this->m_i++; 
}


// == Row cursor

template <typename Scalar>
DenseRowMajor<Scalar>::RowCursor::RowCursor(DenseRowMajor &m)
	: m_matrix(m)
{
}


template <typename Scalar>
void
DenseRowMajor<Scalar>::RowCursor::reset()
{
	this->m_i = 0;
	this->m_j = 0;
}


template <typename Scalar>
Scalar
DenseRowMajor<Scalar>::RowCursor::getElement() const
{
	return this->m_matrix.getValue(this->m_i,this->m_j);
}


template <typename Scalar>
void
DenseRowMajor<Scalar>::RowCursor::setElement(Scalar value) const
{
	this->m_matrix.setValue(this->m_i,this->m_j, value);
}


template <typename Scalar>
void
DenseRowMajor<Scalar>::RowCursor::getNextNonNullRowElement() 
{
	this->m_j++;
}


template <typename Scalar>
bool
DenseRowMajor<Scalar>::RowCursor::hasNextNonNullRowElement() const
{
	return this->m_j < this->m_matrix.columns();
}


template <typename Scalar>
void 
DenseRowMajor<Scalar>::RowCursor::startNextNonNullRow() 
{
	this->m_i = 0;
	this->m_j++;
}


template <typename Scalar>
bool
DenseRowMajor<Scalar>::RowCursor::hasNextNonNullRow() const
{
	return this->m_i < this->m_matrix.rows();
}


// Specify traits

template<typename Scalar>
struct Traits<DenseRowMajor<Scalar> >
{
	static constexpr bool is_writeable() noexcept { return true; }
	static constexpr bool is_resizeable() noexcept { return true; }
	static constexpr bool is_sparse() noexcept { return false; }
};


}	// namespace matrix
}	// namespace mla
#endif
