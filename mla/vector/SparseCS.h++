#ifndef MLA_VECTOR_STORAGE_POLICY_SPARSE_CS_HPP
#define MLA_VECTOR_STORAGE_POLICY_SPARSE_CS_HPP


#include <mla/vector/traits.h++>
#include <mla/LAException.h++>
#include <vector>

#include <mla/VectorCursor.h++>


namespace mla
{
namespace vector
{

/**
SparseCS: a storage policy class intended to implement the interfaces specific 
to the sparse, dynamically resizeable vector 
**/
template<typename Scalar>
class SparseCS
{
public:
	typedef Scalar scalar_type;

	struct Data
	{
		size_t	t_size;		// holds the vector's size

		std::vector< size_t >	column_index;	
		std::vector< Scalar >	values;	
	} data;
	
public:
	SparseCS(size_t size = 0);

	/**
	 * Returns the size of the vector
	 */
	size_t size() const		{ return data.t_size; };

	/**
	 * Sets all values to zero
	 */
	void setZero();

	/**
	 * Returns the value in [index]
	 */
	Scalar getValue(size_t index) const;

	/**
	 * Sets the value at index. If the element doesn't exist, inserts it.
	 * @param index	index of the vector element
	 * @param val	value to be assigned to the element's vector
	 */
	void setValue(size_t index, Scalar val);

	/**
	 * Returns a reference to the element in (row, column)
	 * @param index	index of the vector element
	 */
	Scalar & operator() (size_t index);

	/**
	 * Changes the size of the vector, resetting it in the process.
	 * @param size	new vector size
	 */
	void resize(size_t size);

	/**
	 * Reserve memory for at least size-many non-zero entries
	 * @param	size	number of entries to reserve memory for
	 */
	void reserve(size_t size);

	void push_back(size_t index, Scalar value);


	/**
	 * The cursor that navigates through the non-null elements
	 * of this vector
	 */
	class Cursor
		: public VectorCursor<Scalar, SparseCS>
	{
	public:
		Cursor(SparseCS<Scalar> &m);

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
	Cursor cursor()	{	return	typename SparseCS<Scalar>::Cursor(*this);	}

};



template<typename Scalar>
SparseCS<Scalar>::SparseCS(size_t size)
{
	resize(size);
}


template<typename Scalar>
void 
SparseCS<Scalar>::setZero()
{
	data.values.clear();
	data.column_index.clear();
}


template<typename Scalar>
Scalar 
SparseCS<Scalar>::getValue(size_t index) const
{
	if(index >= data.t_size)
	{
		throw LAException("SparseCS::getValue() index >= data.t_size");
	}

	for(size_t i = 0; i != data.column_index.size(); i++)
	{
		if(data.column_index[i] == index)
			return data.values[i];
		if(data.column_index[i] > index)
			return 0;
	}

	return 0;
}


template<typename Scalar>
void
SparseCS<Scalar>::setValue (size_t index, Scalar value)
{
	if(index >= data.t_size)
	{
		throw LAException("SparseCS::operator() index >= data.t_size");
	}

	for(size_t i = 0; i < data.column_index.size(); i++)
	{
		if(data.column_index[i] == index)
		{
			data.values[i] = value;
			return;
		}

		if(data.column_index[i] > index)
		{
			std::vector<size_t>::iterator ci = data.column_index.begin() + i;
			data.column_index.insert(ci,index);
			typename std::vector<Scalar>::iterator cv = data.values.begin() + i;
			data.values.insert(cv, value);
		}
	}

	data.column_index.push_back(index);
	data.values.push_back( value );
}


template<typename Scalar>
Scalar &
SparseCS<Scalar>::operator() (size_t index)
{
	if(index >= data.t_size)
	{
		throw LAException("SparseCS::operator() index >= data.t_size");
	}

	size_t i = 0;
	for(i = 0; i < data.column_index.size(); i++)
	{
		if(data.column_index[i] == index)
		{
			break;
		}

		if(data.column_index[i] > index)
		{
			std::vector<size_t>::iterator ci = data.column_index.begin() + i;
			data.column_index.insert(ci,index);
			typename std::vector<Scalar>::iterator cv = data.values.begin() + i;
			data.values.insert(cv, (Scalar) 0);
		}
	}

	return data.values[i];
}


template<typename Scalar>
void 
SparseCS<Scalar>::resize(size_t size)
{
	using namespace std;

	data.column_index.clear();
	data.values.clear();
	data.t_size = size;
}


template<typename Scalar>
void 
SparseCS<Scalar>::reserve(size_t size)
{
	using namespace std;

	data.column_index.reserve(size);
	data.values.reserve(size);
}


template<typename Scalar>
void 
SparseCS<Scalar>::push_back(size_t size, Scalar value)
{
	using namespace std;

	data.column_index.push_back(size);
	data.values.push_back(value);
}


/**
 * === The Cursor implementation
 **/

template<typename Scalar>
SparseCS<Scalar>::Cursor::Cursor(SparseCS<Scalar> &m)
	: VectorCursor<Scalar, SparseCS> (m)
{
	this->reset();
}


template<typename Scalar>
void 
SparseCS<Scalar>::Cursor::reset()
{
	if(this->m_vector.data.values.size() == 0)
		this->m_i = this->m_vector.data.column_index.size();
	else
		this->m_i = 0;
}


template<typename Scalar>
Scalar 
SparseCS<Scalar>::Cursor::element() const
{
	return this->m_vector.data.values[this->m_i];
}


template<typename Scalar>
size_t 
SparseCS<Scalar>::Cursor::current() const
{
	return this->m_vector.data.column_index[this->m_i];
}


template<typename Scalar>
bool 
SparseCS<Scalar>::Cursor::at_end() const
{
	return this->m_i >= this->m_vector.data.column_index.size();
}


template<typename Scalar>
void 
SparseCS<Scalar>::Cursor::next()
{
	this->m_i++;
}


// vector traits


template<typename Scalar>
class Traits<SparseCS<Scalar> >
{
	typedef Scalar scalar_type;

	static constexpr bool is_writeable() noexcept { return true; }
	static constexpr bool is_resizeable() noexcept { return true; }
	static constexpr bool is_sparse() noexcept { return true; }
};

}	// namespace mla::vector

}	// namespace mla

#endif
