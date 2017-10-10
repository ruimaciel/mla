#ifndef MLA_MATRIX_STORAGE_POLICY_SPARSE_DOK_HPP
#define MLA_MATRIX_STORAGE_POLICY_SPARSE_DOK_HPP

#include <map>

#include <mla/matrix/traits.h++>
#include <mla/LAException.h++>
#include <mla/vector/SparseCS.h++>


namespace mla
{
namespace matrix
{

/**
SparseDOK: a storage policy class for the Matrix host class.
This class implements the interface for the sparse Dictionary of Keys (DOK) 
matrix format
**/
template<typename Scalar>
class SparseDOK
{
public:
	typedef Scalar scalar_type;

	// type used in row and column vectors
	typedef vector::SparseCS<Scalar>	VectorType;

	struct Data
	{
		size_t	n_rows;		// number of rows
		size_t	n_columns;	// number of columns

		std::map< std::pair<size_t,size_t>, Scalar> 	key_value_map;
		Scalar & getKeyReference(size_t i, size_t j)	{	return this->key_value_map[std::pair<size_t,size_t>(i,j)];	}
	} data;

public:
	SparseDOK(size_t rows = 1, size_t columns = 1);

	/**
	Sets all values to zero
	**/
	void setZero();

	size_t rows() const		{ return data.n_rows; };
	size_t columns() const		{ return data.n_columns; };

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
	* Returns the number of non-zero coefficients
	**/
	size_t nnz() const;


	/**
	 * Permutes two rows
	 *@param row1 
	 *@param row2 
	 **/
	void permuteRows(size_t row1, size_t row2);


	/**
	 * Permutes two columns
	 *@param column1 
	 *@param column2 
	 **/
	void permuteColumns(size_t column1, size_t column2);


	/**
	 * Reorders rows and columns
	 **/
	void symmetricReorder(std::vector<size_t> const & indices);
};



template<typename Scalar>
SparseDOK<Scalar>::SparseDOK(size_t rows, size_t columns)
{
	resize(rows, columns);
}


template<typename Scalar>
void
SparseDOK<Scalar>::setZero()
{
	data.key_value_map.clear();
}


template<typename Scalar>
Scalar
SparseDOK<Scalar>::getValue(size_t row, size_t column) const
{
	using namespace std;

	if(row >= data.n_rows)
	{
		throw LAException("row >= rows()");
	}
	if(column >= data.n_columns)
	{
		throw LAException("column >= columns()");
	}

	typename map< std::pair<size_t,size_t>, Scalar>::const_iterator	i;

	std::pair<size_t, size_t> key(row,column);
	i = data.key_value_map.find(key);
	if(i == data.key_value_map.end())
		return (Scalar)0;

	return i->second;
}



template<typename Scalar>
void
SparseDOK<Scalar>::setValue(size_t row, size_t column, Scalar value)
{
	if(row >= data.n_rows)
	{
		throw LAException("row >= rows()");
	}
	if(column >= data.n_columns)
	{
		throw LAException("column >= columns()");
	}


	this->data.getKeyReference(row,column) = value;
}



template<typename Scalar>
Scalar & 
SparseDOK<Scalar>::operator() (size_t row, size_t column)
{
	if(row >= data.n_rows)
	{
		throw LAException("row >= rows()");
	}
	if(column >= data.n_columns)
	{
		throw LAException("column >= columns()");
	}

	return this->data.getKeyReference(row,column);
}


template<typename Scalar>
void
SparseDOK<Scalar>::resize(size_t rows, size_t columns)
{
	using namespace std;

	data.key_value_map.clear();

	data.n_rows = rows;
	data.n_columns = columns;
}


template<typename Scalar>
void
SparseDOK<Scalar>::setEye()
{
	data.key_value_map.clear();
	
	auto n_elements = rows()< columns() ? rows(): columns();

	for(size_t k = 0; k < n_elements; k++)
	{
		this->setValue(k,k,(Scalar)1);
	}
}


template<typename Scalar>
size_t 
SparseDOK<Scalar>::nnz() const
{
	return this->data.key_value_map.size();
}


template<typename Scalar>
void 
SparseDOK<Scalar>::permuteRows(size_t row1, size_t row2)
{
	typename std::map< std::pair<size_t, size_t>, Scalar> row1_elements, row2_elements;

	// Extract both rows 
	for(auto it = data.key_value_map.begin(); it != data.key_value_map.end(); it++)
	{
		auto current_key = it->first;
		if( current_key.first == row1)
		{
			current_key.first = row2;
			row2_elements[current_key] = it->second;
			data.key_value_map.erase(it);
		}
		else if( current_key.first == row2)
		{
			current_key.first = row1;
			row1_elements[current_key] = it->second;
			data.key_value_map.erase(it);
		}
	}

	data.key_value_map.insert(row1_elements.begin(), row1_elements.end());
	row1_elements.clear();	// free memory
	data.key_value_map.insert(row2_elements.begin(), row2_elements.end());
}


template<typename Scalar>
void 
SparseDOK<Scalar>::permuteColumns(size_t column1, size_t column2)
{
	typename std::map< std::pair<size_t, size_t>, Scalar> column1_elements, column2_elements;

	// Extract both rows 
	for(auto it = data.key_value_map.begin(); it != data.key_value_map.end(); it++)
	{
		auto current_key = it->first;
		if( current_key.second == column1)
		{
			current_key.second = column2;
			column2_elements[current_key] = it->second;
			data.key_value_map.erase(it);
		}
		else if( current_key.second == column2)
		{
			current_key.second = column1;
			column1_elements[current_key] = it->second;
			data.key_value_map.erase(it);
		}
	}

	data.key_value_map.insert(column1_elements.begin(), column1_elements.end());
	column1_elements.clear();	// free memory
	data.key_value_map.insert(column2_elements.begin(), column2_elements.end());
}



template<typename Scalar>
void 
SparseDOK<Scalar>::symmetricReorder(std::vector<size_t> const & indices)
{
	// check if indices are compatible with the matrix
	if(this->columns() != indices.size())
	{
		throw LAException("indices must match matrix size");
	}

	typename std::map< std::pair<size_t,size_t>, Scalar> new_kvm;
	for( typename std::map< std::pair<size_t,size_t>, Scalar>::iterator it = data.key_value_map.begin(); it != data.key_value_map.end();)
	{
		std::pair<size_t, size_t> const & old_key = it->first;

		// prepare new key
		std::pair<size_t, size_t> new_key;
		new_key.first = indices[old_key.first];
		new_key.second = indices[old_key.second];

		new_kvm[new_key] = it->second;

		// 
		data.key_value_map.erase(it++);
	}

	// updates the key
	data.key_value_map = new_kvm;
}


// Specify traits

template<typename Scalar>
struct Traits<SparseDOK<Scalar> >
{
	static constexpr bool is_writeable() noexcept { return true; }
	static constexpr bool is_resizeable() noexcept { return true; }
	static constexpr bool is_sparse() noexcept { return true; }
};


}	// namespace matrix
}	// namespace mla



#endif
