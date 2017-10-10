#ifndef MLA_VECTOR_CURSOR_HPP
#define MLA_VECTOR_CURSOR_HPP

namespace mla
{



/**
 * Base class that defines the interface used by all matrix cursor libraries
 */
template<typename Scalar, template<typename> class VectorStoragePolicy >
class VectorCursor
{
protected:
	VectorStoragePolicy<Scalar> &m_vector;

	size_t m_i;	// current row

public:
	VectorCursor(VectorStoragePolicy<Scalar> &m);

	/**
	 * resets the cursor to the first non-null term of the matrix
	 */
	virtual void reset() = 0;

	/**
	 * returns the element
	 **/
	virtual Scalar element() const = 0;

	/**
	 * Returns the row currently being pointed out by the cursor
	 */
	virtual size_t current() const = 0;

	/**
	 * Check if cursor is pointing to the element following the last one of the current column
	 */
	virtual bool at_end() const = 0;

	virtual void next() = 0;
};



template<typename Scalar, template<typename> class VectorStoragePolicy >
VectorCursor<Scalar, VectorStoragePolicy>::VectorCursor( VectorStoragePolicy<Scalar> &m)
	: m_vector(m)
{
}


}	// namespace mla

#endif
