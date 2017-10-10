#ifndef MLA_MATRIX_STORAGE_POLICY_SPARSE_CCS_HPP
#define MLA_MATRIX_STORAGE_POLICY_SPARSE_CCS_HPP

#include <vector>

#include <mla/matrix/traits.h++>
#include <mla/LAException.h++>
#include <mla/vector/SparseCS.h++>



namespace mla
{
namespace matrix
{

/**
SparseCCS: a storage policy class for the Matrix host class.
This class implements the interface for the sparse Compressed Row Storage (CRS) sparse matrix format
http://netlib.org/linalg/html_templates/node91.html#SECTION00931100000000000000
**/
template<typename Scalar>
class SparseCCS
{
public:
	typedef Scalar scalar_type;

	// type used in row and column vectors
	typedef vector::SparseCS<Scalar>	VectorType;

	struct Data
	{
		size_t	n_rows;	// number of columns

		std::vector<Scalar> values;
		std::vector<size_t> row_index;
		std::vector<size_t> column_pointer;
	} data;


public:
	SparseCCS(size_t rows = 0, size_t columns = 0);

	/**
	Sets all values to zero
	**/
	void setZero();

	size_t rows() const		{ return data.n_rows; };
	size_t columns() const		{ return data.column_pointer.size()-1; };

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
};



template<typename Scalar>
SparseCCS<Scalar>::SparseCCS(size_t rows, size_t columns)
{
	resize(rows, columns);
}



template<typename Scalar>
void
SparseCCS<Scalar>::setZero()
{
	//TODO properly clear the sparse matrix instead of just zeroing the existing values
	for(typename std::vector<Scalar>::iterator i = data.values.begin(); i != data.values.end(); i++)
	{
		*i = 0;
	}
}


template<typename Scalar>
Scalar
SparseCCS<Scalar>::getValue(size_t row, size_t column) const
{
	using namespace std;

	if(row >= this->rows() )
	{
		throw LAException("row >= rows()");
	}
	if(column >= this->columns() )
	{
		throw LAException("column >= columns()");
	}

	size_t i = data.column_pointer[column];

	if(row == data.row_index[i])
		return data.values[i];

	if(row < data.row_index[i])
		return 0;

	for(i = data.column_pointer[column]; i < data.column_pointer[column+1]; i++)
	{
		if(data.row_index[i] < row)
			continue;

		if(data.row_index[i] == row)	
			return data.values[i];

		break;
	}
	return 0;
}



template<typename Scalar>
void
SparseCCS<Scalar>::setValue(size_t row, size_t column, Scalar value)
{
	if(row >= (size_t)this->rows() )
	{
		throw LAException("row >= rows()");
	}
	if(column >= (size_t)this->columns() )
	{
		throw LAException("column >= columns()");
	}


	if(row == (size_t)data.row_index[data.column_pointer[column]])
	{
		data.values[ data.column_pointer[column]] = value;
		return;
	}

	// the first element in this column has a greater row index than the one being referenced
	if(row < (size_t)data.row_index[data.column_pointer[column]])
	{
		std::vector<size_t>::iterator ro = data.row_index.begin();
		advance(ro,data.column_pointer[column]);
		data.row_index.insert(ro,row);

		typename std::vector<Scalar>::iterator val = data.values.begin();
		advance(val, data.column_pointer[column]);
		data.values.insert(val, 0);

		for(size_t i = column+1; i < data.column_pointer.size(); i++)
		{
			data.column_pointer[i]++;
		}
		
		data.values[ data.column_pointer[column]] = value;
		return;
	}

	// the first element in this row has an inferior column index.  Search for the first superior instance, insert a new element and reference it
	size_t j;
	for(j = data.column_pointer[column]; j < data.column_pointer[column+1]; j++)
	{
		if((size_t)data.row_index[j] < row)
			continue;

		if((size_t)data.row_index[j] == row)	
		{
			data.values[j] = value;
			return;
		}

		break;
	}

	std::vector<size_t>::iterator ro = data.row_index.begin();
	advance(ro,j);
	data.row_index.insert(ro,row);

	typename std::vector<Scalar>::iterator val = data.values.begin();
	advance(val, j);
	data.values.insert(val, 9);

	for(size_t i = column+1; i < data.column_pointer.size(); i++)
	{
		data.column_pointer[i]++;
	}

	data.values[ j] = value;
}



template<typename Scalar>
Scalar &
SparseCCS<Scalar>::operator() (size_t row, size_t column)
{
	//TODO finish this
	using namespace std;

	if(row >= (size_t)this->rows() )
	{
		throw LAException("row >= rows()");
	}
	if(column >= (size_t)this->columns() )
	{
		throw LAException("column >= columns()");
	}


	if(row == (size_t)data.row_index[data.column_pointer[column]])
		return data.values[ data.column_pointer[column]];

	// the first element in this column has a greater row index than the one being referenced
	if(row < (size_t)data.row_index[data.column_pointer[column]])
	{
		std::vector<size_t>::iterator ro = data.row_index.begin();
		advance(ro,data.column_pointer[column]);
		data.row_index.insert(ro,row);

		typename std::vector<Scalar>::iterator val = data.values.begin();
		advance(val, data.column_pointer[column]);
		data.values.insert(val, 0);

		for(size_t i = column+1; i < data.column_pointer.size(); i++)
		{
			data.column_pointer[i]++;
		}
		
		return data.values[ data.column_pointer[column]];
	}

	// the first element in this row has an inferior column index.  Search for the first superior instance, insert a new element and reference it
	size_t j;
	for(j = data.column_pointer[column]; j < data.column_pointer[column+1]; j++)
	{
		if((size_t)data.row_index[j] < row)
			continue;

		if((size_t)data.row_index[j] == row)	
			return data.values[j];

		break;
	}

	std::vector<size_t>::iterator ro = data.row_index.begin();
	advance(ro,j);
	data.row_index.insert(ro,row);

	typename std::vector<Scalar>::iterator val = data.values.begin();
	advance(val, j);
	data.values.insert(val, 9);

	for(size_t i = column+1; i < data.column_pointer.size(); i++)
	{
		data.column_pointer[i]++;
	}

	return data.values[ j];
}


template<typename Scalar>
void
SparseCCS<Scalar>::resize(size_t rows, size_t columns)
{
	using namespace std;

	//TODO find a better way to initialize a CRS matrix than setting zeros at the diagonal

	data.values.resize(columns);
	data.row_index.resize(columns);
	data.column_pointer.resize(columns+1);	// to account to the after-last element guard

	data.n_rows = rows;

	for(size_t i = 0; i < columns; i++)
	{
		data.row_index[i] = i;
		data.values[i] = (Scalar)0;
	}

	for(size_t j = 0; j <= columns; j++)
		data.column_pointer[j] = j;
}


template<typename Scalar>
void
SparseCCS<Scalar>::setEye()
{
	auto n_elements = rows() < columns()? rows(): columns();

	data.values.resize(n_elements,(Scalar)1);
	data.row_index.resize(n_elements);

	size_t k;
	for(k = 0; k < n_elements; k++)
	{
		data.values[k] = (Scalar)1;
		data.row_index[k] = k;
		data.column_pointer[k] = k;
	}

	for(size_t j = k; j < columns(); j++)
	{
		data.column_pointer[j] = k;
	}
}


// Specify traits

template<typename Scalar>
struct Traits<SparseCCS<Scalar> >
{
	static constexpr bool is_writeable() noexcept { return true; }
	static constexpr bool is_resizeable() noexcept { return true; }
	static constexpr bool is_sparse() noexcept { return true; }
};



}	// namespace matrix
}	// namespace mla
#endif
