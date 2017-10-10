#ifndef MLA_MATRIX_STORAGE_POLICY_SPARSE_CRS_HPP
#define MLA_MATRIX_STORAGE_POLICY_SPARSE_CRS_HPP

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
SparseCRS: a storage policy class for the Matrix host class.
This class implements the interface for the sparse Compressed Row Storage (CRS) sparse matrix format
http://netlib.org/linalg/html_templates/node91.html#SECTION00931100000000000000
**/
template<typename Scalar>
class SparseCRS
{
public:
	typedef Scalar scalar_type;

	// type used in row and column vectors
	typedef vector::SparseCS<Scalar>	VectorType;

	struct Data 
	{
		size_t	n_columns;	// number of columns

		std::vector<Scalar> values;
		std::vector<size_t> column_index;
		std::vector<size_t> row_pointer;
	} data;


public:
	SparseCRS(size_t rows = 0, size_t columns = 0);

	/**
	Sets all values to zero
	**/
	void setZero();

	size_t rows() const		{ return data.row_pointer.size()-1; };
	size_t columns() const		{ return data.n_columns; };

	/*
	Returns the value in [row,column] 
	*/
	Scalar getValue(size_t row, size_t column) const;

	/**
	 * Returns the row vector
	 **/
	VectorType getRow(size_t row) const;

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
		: public MatrixCursor<Scalar, SparseCRS>
	{
	public:
		Cursor(SparseCRS &m);

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
	Cursor cursor()	{	return	typename SparseCRS<Scalar>::Cursor(*this);	}


	/**
	 * A cursor design pattern for row-wise iterations
	 **/
	class RowCursor 
		: public matrix::RowCursor<Scalar>
	{
	protected:
		SparseCRS &m_matrix;
		size_t m_ci, m_rp;

	public:
		RowCursor(SparseCRS &matrix);

		void reset() override;
		Scalar getElement() const override;
		void setElement(Scalar value) const override;

		size_t getCurrentRow() const override;

		/**
		 * Returns the column currently being pointed out by the cursor
		 */
		size_t getCurrentColumn() const override;

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
	RowCursor	getRowCursor()	{	return	typename SparseCRS<Scalar>::RowCursor(*this);	}
};



template<typename Scalar>
SparseCRS<Scalar>::SparseCRS(size_t rows, size_t columns)
{
	resize(rows, columns);
}


template<typename Scalar>
void
SparseCRS<Scalar>::setZero()
{
	data.values.resize(1);
	data.values[0] = (Scalar)0;

	data.column_index.resize(1);
	data.column_index[0] = 0;

	data.row_pointer.resize(this->rows()+1);
	data.row_pointer[0] = 0;

	for(size_t k = 1; k < data.row_pointer.size(); k++)
	{
		data.row_pointer[k] = 1;
	}
}


template<typename Scalar>
Scalar
SparseCRS<Scalar>::getValue(size_t row, size_t column) const
{
	if(row >= this->rows())
	{
		throw LAException("row >= rows()");
	}
	if(column >= this->columns())
	{
		throw LAException("column >= columns()");
	}

	for(size_t i = data.row_pointer[row]; i < data.row_pointer[row+1]; i++)
	{
		if(data.column_index[i] == column)	
		{
			return data.values[i];
		}
		else if(data.column_index[i] > column)
			break;
	}
	return (Scalar)0;
}


template<typename Scalar>
typename SparseCRS<Scalar>::VectorType 
SparseCRS<Scalar>::getRow(size_t row) const
{
	if(row >= this->rows() )
	{
		throw LAException("row >= rows()");
	}

	using VectorType = typename SparseCRS<Scalar>::VectorType;

	VectorType row_vector(this->columns());

	// get elements 
	size_t start_index = this->data.row_pointer[row];
	size_t end_index = this->data.row_pointer[row+1];

	row_vector.reserve(end_index-start_index);

	for(size_t j = start_index; j < end_index; j++)
	{
		size_t const &column_index = this->data.column_index[j];
		Scalar const &value = this->data.values[j];
		row_vector.push_back(column_index, value);
	}

	return row_vector;
}


template<typename Scalar>
void
SparseCRS<Scalar>::setValue(size_t row, size_t column, Scalar value)
{
	//TODO finish this

	if(row >= this->rows())
	{
		throw LAException("row >= rows()");
	}
	if(column >= this->columns())
	{
		throw LAException("column >= columns()");
	}


	if(column == data.column_index[data.row_pointer[row]])
	{
		data.values[ data.row_pointer[row]] = value;
	}

	// the first element in this row has a greater column index than the one being referenced
	if(data.column_index[data.row_pointer[row]] > column)
	{
		std::vector<size_t>::iterator col = data.column_index.begin();
		advance(col,data.row_pointer[row]);
		data.column_index.insert(col,column);

		typename std::vector<Scalar>::iterator val = data.values.begin();
		advance(val, data.row_pointer[row]);
		data.values.insert(val, 0);

		for(size_t i = row+1; i < data.row_pointer.size(); i++)
		{
			data.row_pointer[i]++;
		}
		
		data.values[ data.row_pointer[row]] = value;
		return;
	}

	// the first element in this row has an inferior column index.  Search for the first superior instance, insert a new element and reference it
	size_t j;
	for(j = data.row_pointer[row]; j < data.row_pointer[row+1]; j++)
	{
		if(data.column_index[j] < column)
			continue;

		if(data.column_index[j] == column)	
		{
			data.values[j] = value;
			return;
		}

		break;
	}

	std::vector<size_t>::iterator col = data.column_index.begin();
	advance(col,j);
	data.column_index.insert(col,column);

	typename std::vector<Scalar>::iterator val = data.values.begin();
	advance(val, j);
	data.values.insert(val, 0);

	for(size_t i = row+1; i < data.row_pointer.size(); i++)
	{
		data.row_pointer[i]++;
	}

	data.values[ j] = value;

}



template<typename Scalar>
Scalar &
SparseCRS<Scalar>::operator() (size_t row, size_t column)
{
	//TODO finish this

	if(row >= this->rows())
	{
		throw LAException("row >= rows()");
	}
	if(column >= this->columns())
	{
		throw LAException("column >= columns()");
	}


	// needs to search for the right element in the current row
	size_t idx = data.row_pointer[row];
	for(; idx < data.row_pointer[row+1]; idx++)
	{
		size_t j = data.column_index[idx];
		if(j == column)
		{
			// NNZ element was found.  just return a reference to it
			return data.values[idx];
		}
		else if (j > column)
		{
			// NNZ element doesn't exist
			break;
		}
	}

	// element wasn't added.  Let's add a NNZ element
	if(idx == data.values.size())
	{
		data.column_index.push_back(column);
		data.values.push_back( (Scalar)0);
	}
	else
	{
		std::vector<size_t>::iterator col_j = data.column_index.begin();
		std::advance(col_j, idx);
		data.column_index.insert(col_j,column);

		typename std::vector<Scalar>::iterator nnz_j = data.values.begin();
		std::advance(nnz_j, idx);
		data.values.insert(nnz_j, (Scalar)0);
	}

	// let's update the column index
	std::vector<size_t>::iterator row_ptr_iter = data.row_pointer.begin();
	std::advance(row_ptr_iter, row+1);
	for(; row_ptr_iter != data.row_pointer.end(); row_ptr_iter++)
	{
		(*row_ptr_iter)++;
	}

	return data.values[ idx];

}


template<typename Scalar>
void
SparseCRS<Scalar>::resize(size_t rows, size_t columns)
{
	data.values.resize(1);
	data.column_index.resize(1);
	data.row_pointer.resize(rows+1);

	data.values[0] = (Scalar)0;
	data.column_index[0] = 0;
	data.row_pointer[0] = 0;

	for(size_t i = 1; i < rows; i++)
	{
		data.row_pointer[i] = 1;
	}

	data.n_columns = columns;
}


template<typename Scalar>
void
SparseCRS<Scalar>::setEye()
{
	auto n_elements = rows() < columns()? rows(): columns();

	data.values.resize(n_elements,(Scalar)1);
	data.column_index.resize(n_elements);

	size_t k;
	for(k = 0; k < n_elements; k++)
	{
		data.values[k] = (Scalar)1;
		data.column_index[k] = k;
		data.row_pointer[k] = k;
	}

	for(size_t i = k; i < data.row_pointer.size(); i++)
	{
		data.row_pointer[i] = k;
	}
}


/**
 * Start of cursor definition
 **/

template<typename Scalar>
SparseCRS<Scalar>::Cursor::Cursor(SparseCRS &m)
	: MatrixCursor<Scalar, SparseCRS> (m)
{
	this->reset();
}


template<typename Scalar>
void
SparseCRS<Scalar>::Cursor::start_current_column()
{
	this->m_j = 0;
}


template<typename Scalar>
void
SparseCRS<Scalar>::Cursor::start_current_row()
{
	this->m_i = 0;
}


template<typename Scalar>
void
SparseCRS<Scalar>::Cursor::start_current_column_nn()
{
	auto & data = this->m_matrix.data;
	size_t i = 0;
	for(; i < data.column_index.size(); i++)
	{
		if( data.column_index[i] == this->m_j)
			break;
	}

	if(i >= data.column_index.size())
	{
		// point to EOR
		this->m_i = this->m_matrix.rows();
		return;
	}

	// i points to column between 
	this->m_i = i;
}


template<typename Scalar>
void
SparseCRS<Scalar>::Cursor::start_current_row_nn()
{
	size_t col_index = this->m_matrix.data.row_pointer[this->m_i];
	this->m_j = this->m_matrix.data.column_index[col_index];
}


template<typename Scalar>
void
SparseCRS<Scalar>::Cursor::reset()
{
	this->m_i = 0;
	this->m_j = 0;
}


template<typename Scalar>
Scalar
SparseCRS<Scalar>::Cursor::element()
{
	return this->m_matrix.getValue( this->current_row(), this->current_column() );
}


template<typename Scalar>
size_t
SparseCRS<Scalar>::Cursor::current_row() const
{
	return this->m_i;
}


template<typename Scalar>
size_t
SparseCRS<Scalar>::Cursor::current_column() const
{
	return this->m_j;
}



template<typename Scalar>
bool
SparseCRS<Scalar>::Cursor::at_beginning_of_column() const
{
	return this->m_j == 0;
}

	
template<typename Scalar>
bool
SparseCRS<Scalar>::Cursor::at_beginning_of_row() const
{
	return this->m_j == 0;
}

template<typename Scalar>
bool
SparseCRS<Scalar>::Cursor::at_end_of_current_column() const
{
	return !(this->m_j < this->m_matrix.columns() );
}

	
template<typename Scalar>
bool
SparseCRS<Scalar>::Cursor::at_end_of_current_row() const
{
	return !(this->m_i < this->m_matrix.rows() );
}


template<typename Scalar>
bool
SparseCRS<Scalar>::Cursor::at_end_of_columns() const
{
	return !(this->m_j < this->m_matrix.columns() );
}

	
template<typename Scalar>
bool
SparseCRS<Scalar>::Cursor::at_end_of_rows() const
{
	return !(this->m_i < this->m_matrix.rows() );
}



template<typename Scalar>
void
SparseCRS<Scalar>::Cursor::start_next_column_nn()
{
	this->m_j++;
	this->start_current_column_nn();
}


template<typename Scalar>
void
SparseCRS<Scalar>::Cursor::start_next_row_nn()
{
	this->m_i++;
	if(this->m_i < this->m_matrix.rows())
	{
		this->m_j = this->m_matrix.data.row_pointer[this->m_i];
	}
	else
	{
		this->m_j = this->m_matrix.columns();
	}
	
}


template<typename Scalar>
void
SparseCRS<Scalar>::Cursor::increment_column()
{
	this->m_j++; 
}


template<typename Scalar>
void
SparseCRS<Scalar>::Cursor::increment_row()
{
	this->m_i++; 
}


// Row Cursor implementaton


template<typename Scalar>
SparseCRS<Scalar>::RowCursor::RowCursor(SparseCRS &matrix)
	: m_matrix(matrix)
{
	this->reset();	
}


template<typename Scalar>
void 
SparseCRS<Scalar>::RowCursor::RowCursor::reset()
{
	this->m_ci = 0;
	this->m_rp = 0;
}


template<typename Scalar>
Scalar 
SparseCRS<Scalar>::RowCursor::RowCursor::getElement() const
{
	return this->m_matrix.data.values[this->m_ci];
}


template<typename Scalar>
void 
SparseCRS<Scalar>::RowCursor::RowCursor::setElement(Scalar value) const
{
	this->m_matrix.data.values[this->m_ci] = value;
}


template<typename Scalar>
size_t 
SparseCRS<Scalar>::RowCursor::getCurrentRow() const
{ 
	return this->m_matrix.data.row_pointer[this->m_rp];
}

template<typename Scalar>
size_t 
SparseCRS<Scalar>::RowCursor::getCurrentColumn() const
{ 
	return this->m_matrix.data.column_index[this->m_ci];
}


template<typename Scalar>
void
SparseCRS<Scalar>::RowCursor::getNextNonNullRowElement()
{
	this->m_ci++;
}


template<typename Scalar>
bool
SparseCRS<Scalar>::RowCursor::hasNextNonNullRowElement() const
{
	return this->m_ci < this->m_matrix.data.row_pointer[this->m_rp];
}


template<typename Scalar>
void
SparseCRS<Scalar>::RowCursor::startNextNonNullRow()
{
	this->m_rp++;
	this->m_ci = this->m_matrix.data.row_pointer[this->m_rp];
}


template<typename Scalar>
bool
SparseCRS<Scalar>::RowCursor::hasNextNonNullRow()const
{
	return this->m_ci < this->m_matrix.data.row_pointer.size();
}


// Specify traits

template<typename Scalar>
struct Traits<SparseCRS<Scalar> >
{
	static constexpr bool is_writeable() noexcept { return true; }
	static constexpr bool is_resizeable() noexcept { return true; }
	static constexpr bool is_sparse() noexcept { return true; }
};


}	// namespace matrix
}	// namespace mla

#endif
