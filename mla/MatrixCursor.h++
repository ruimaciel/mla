#ifndef MLA_MATRIX_CURSOR_HPP
#define MLA_MATRIX_CURSOR_HPP

namespace mla
{



/**
 * Base class that defines the interface used by all matrix cursor libraries
 */
template<typename Scalar, template<typename> class MatrixStoragePolicy >
class MatrixCursor
{
protected:
	MatrixStoragePolicy<Scalar> &m_matrix;

	size_t m_i;	// current row
	size_t m_j;	// current column

public:
	MatrixCursor(MatrixStoragePolicy<Scalar> &m);

	/**
	 * resets the cursor to the first non-null term of the matrix
	 */
	virtual void reset() = 0;

	/**
	 * returns the element
	 **/
	virtual Scalar element() = 0;

	/**
	 * Returns the row currently being pointed out by the cursor
	 */
	virtual size_t current_row() const = 0;

	/**
	 * Returns the column currently being pointed out by the cursor
	 */
	virtual size_t current_column() const = 0;

	/**
	 * Check if cursor is pointing to the first element of the current column
	 */
	virtual bool at_beginning_of_column() const = 0;
	virtual bool at_beginning_of_row() const = 0;

	/**
	 * Check if cursor is pointing to the element following the last one of the current column
	 */
	virtual bool at_end_of_current_column() const = 0;
	virtual bool at_end_of_current_row() const = 0;

	/**
	 * Check if cursor is pointing to an element pointing to the last column
	 */
	virtual bool at_end_of_columns() const = 0;
	virtual bool at_end_of_rows() const = 0;

	/**
	 * Moves cursor to the first element of the current column
	 */
	virtual void start_current_column() = 0;
	virtual void start_current_row() = 0;

	/**
	 * Moves cursor to the first non-null element of the current column
	 */
	virtual void start_current_column_nn() = 0;
	virtual void start_current_row_nn() = 0;

	/**
	 * Moves cursor to the first element of the subsequent column
	 */
	virtual void start_next_column();
	virtual void start_next_row();

	/**
	 * Moves cursor to the first non-null element of the subsequent column
	 */
	virtual void start_next_column_nn() = 0;
	virtual void start_next_row_nn() = 0;

	/**
	 * Moves cursor to next non-null element in the same row
	 */
	virtual void increment_column() = 0;
	virtual void increment_row() = 0;

};



template<typename Scalar, template<typename> class MatrixStoragePolicy >
MatrixCursor<Scalar, MatrixStoragePolicy>::MatrixCursor( MatrixStoragePolicy<Scalar> &m)
	: m_matrix(m)
{
}


template<typename Scalar, template<typename> class MatrixStoragePolicy >
void 
MatrixCursor<Scalar, MatrixStoragePolicy>::start_next_column()
{
	increment_column();
	start_current_row();
}


template<typename Scalar, template<typename> class MatrixStoragePolicy >
void 
MatrixCursor<Scalar, MatrixStoragePolicy>::start_next_row()
{
	increment_row();
	start_current_column();
}



}	// namespace mla

#endif
