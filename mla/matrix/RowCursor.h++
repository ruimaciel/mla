#ifndef MLA_MATRIX_ROW_CURSOR_HPP
#define MLA_MATRIX_ROW_CURSOR_HPP

#include <cstddef>


namespace mla
{

namespace matrix
{



/**
 * Base class that defines the interface used by all matrix cursor libraries
 */
template<typename Scalar>
class RowCursor
{
public:
	/**
	 * resets the cursor to the first non-null term of the matrix
	 */
	virtual void reset() = 0;

	/**
	 * returns the element
	 **/
	virtual Scalar getElement() const = 0;

	virtual void setElement(Scalar value) const = 0;

	/**
	 * Returns the row currently being pointed out by the cursor
	 */
	virtual size_t getCurrentRow() const = 0;

	/**
	 * Returns the column currently being pointed out by the cursor
	 */
	virtual size_t getCurrentColumn() const = 0;

	/**
	 * Moves cursor to next non-null element in the same row
	 */
	virtual void getNextNonNullRowElement() = 0;
	virtual bool hasNextNonNullRowElement() const = 0;

	virtual void startNextNonNullRow() = 0;
	virtual bool hasNextNonNullRow() const = 0;

};



}	// matrix
}	// namespace mla

#endif
