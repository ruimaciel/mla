#ifndef MLA_VECTOR_STORAGE_POLICY_DENSE_HPP
#define MLA_VECTOR_STORAGE_POLICY_DENSE_HPP


#include <mla/vector/traits.h++>
#include <mla/LAException.h++>
#include <vector>

#include <mla/VectorCursor.h++>


namespace mla
{
namespace vector
{

/**
Dense: a storage policy class intended to implement the interfaces specific to 
the dense, dynamically resizeable vector to be used in the Vector host class.
**/
template<typename Scalar>
class Dense
{
public:
	typedef Scalar scalar_type;

	size_t	t_size;		// number of size

	std::vector< Scalar >	data;	
	
public:
	Dense(size_t size = 0);

	/**
	 * Returns the size of the vector
	 */
	size_t size() const		{ return t_size; };

	/**
	 * Sets all values to zero
	 */
	void setZero();

	/**
	 * Returns the value in [index]
	 */
	Scalar getValue(size_t index) const;

	/**
	 * Sets the value at index.
	 * @param index	index of the vector element
	 * @param val	value to be assigned to the element's vector
	 */
	void setValue(size_t index, Scalar val);

	Scalar & operator[] (size_t index);

	/**
	 * Changes the size of the vector, resetting it in the process.
	 * @param size	new vector size
	 */
	void resize(size_t size);

	/**
	 * The cursor that navigates through the non-null elements
	 * of this vector
	 */
	class Cursor
		: public VectorCursor<Scalar, Dense>
	{
	public:
		Cursor(Dense<Scalar> &m);

		/**
		 * resets the cursor to the first non-null term of the matrix
		 */
		void reset() override final;

		/**
		 * returns the element that the cursor is currently pointing at
		 */
		Scalar element() const override final;

		/**
		 * Returns the row currently being pointed out by the cursor
		 */
		size_t current() const override;

		/**
		 * Check if cursor is pointing to the element following the last one of the current column
		 */
		bool at_end() const override final;

		/**
		 * Moves the cursor to the next non-null element in the vector
		 */
		void next() override final;
	};


	/**
	 * returns a Dense<Scalar>::Cursor object associated with this Dense<Scalar> object
	 */
	Cursor cursor()	{	return	typename Dense<Scalar>::Cursor(*this);	}

};



template<typename Scalar>
Dense<Scalar>::Dense(size_t size)
{
	resize(size);
}


template<typename Scalar>
void 
Dense<Scalar>::setZero()
{
	for(typename std::vector<Scalar>::iterator i = data.begin(); i != data.end(); i++)
	{
		*i = 0;
	}
}


template<typename Scalar>
Scalar 
Dense<Scalar>::getValue(size_t index) const
{
	if(index >= t_size)
	{
		throw LAException("Dense::getValue() index >= data.t_size");
	}

	return this->data[index];
}



template<typename Scalar>
void
Dense<Scalar>::setValue(size_t index, Scalar val)
{
	if(index >= t_size)
	{
		throw LAException("Dense::setValue() index >= data.t_size");
	}

	this->data[index] = val;
}


template<typename Scalar>
Scalar &
Dense<Scalar>::operator[] (size_t index)
{
	if(index >= this->size() )
	{
		throw LAException("Dense::operator[] index >= data.t_size");
	}

	return this->data[index];
}


template<typename Scalar>
void 
Dense<Scalar>::resize(size_t size)
{
	data.resize(size);
	t_size = size;
}


/**
 * === The Cursor implementation
 **/

template<typename Scalar>
Dense<Scalar>::Cursor::Cursor(Dense<Scalar> &m)
	: VectorCursor<Scalar, Dense> (m)
{
	this->reset();
}

template<typename Scalar>
void 
Dense<Scalar>::Cursor::reset()
{
	this->m_i = 0;
}


template<typename Scalar>
Scalar 
Dense<Scalar>::Cursor::element() const
{
	return this->m_vector.data[this->m_i];
}


template<typename Scalar>
size_t 
Dense<Scalar>::Cursor::current() const
{
	return this->m_i;
}


template<typename Scalar>
bool 
Dense<Scalar>::Cursor::at_end() const
{
	return !(this->m_i < this->m_vector.data.size());
}


template<typename Scalar>
void 
Dense<Scalar>::Cursor::next()
{
	this->m_i++;
}



// vector traits


template<typename Scalar>
class Traits<Dense<Scalar> >
{
	typedef Scalar scalar_type;

	static constexpr bool is_writeable() noexcept { return true; }
	static constexpr bool is_resizeable() noexcept { return true; }
	static constexpr bool is_sparse() noexcept { return false; }
};

}	// namespace mla::vector

}	// namespace mla

#endif
