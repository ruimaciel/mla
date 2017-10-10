#ifndef MLA_MATRIX_STORAGE_POLICY_SPARSE_COO_HPP
#define MLA_MATRIX_STORAGE_POLICY_SPARSE_COO_HPP

#include <list>
#include <algorithm>

#include <mla/matrix/traits.h++>
#include <mla/LAException.h++>
#include <mla/MatrixCursor.h++>
#include <mla/vector/SparseCS.h++>


namespace mla
{
namespace matrix
{

/**
SparseCOO: a storage policy class for the Matrix host class.
This class implements the interface for the sparse coordinate list (COO) 
matrix format
**/
template<typename Scalar>
class SparseCOO
{
public:
	typedef Scalar scalar_type;

	// type used in row and column vectors
	typedef vector::SparseCS<Scalar>	VectorType;

	struct Data
	{
		size_t	n_rows;		// number of rows
		size_t	n_columns;	// number of columns

		struct Coordinate {
			size_t i;
			size_t j;
			Scalar value;
		};

		std::list<Coordinate> coordinate_list;

		Scalar & getKeyReference(size_t const i, size_t const j)	
		{	
			auto predicate = [&](typename Data::Coordinate const &c) { return c.i == i && c.j == j; };
			auto it = std::find_if(coordinate_list.begin(), coordinate_list.end(), predicate);
			if( it == coordinate_list.end())
			{
				Coordinate c{i,j,(Scalar)0};
				coordinate_list.push_back(c);
				return coordinate_list.back().value;
			}
			else
			{
				return it->value;
			}
		}
	} data;

public:
	SparseCOO(size_t rows = 1, size_t columns = 1);

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


	/**
	 * returns an object with the bottom part of the current matrix
	 **/
	SparseCOO<Scalar> splitOutBottom(std::size_t row_index);


	/**
	 * returns an object with the right part of the current matrix
	 **/
	SparseCOO<Scalar> splitOutRight(std::size_t column_index);
};



template<typename Scalar>
SparseCOO<Scalar>::SparseCOO(size_t rows, size_t columns)
{
	resize(rows, columns);
}


template<typename Scalar>
void
SparseCOO<Scalar>::setZero()
{
	using namespace std;
	data.coordinate_list.clear();
}


template<typename Scalar>
Scalar
SparseCOO<Scalar>::getValue(size_t row, size_t column) const
{
	using namespace std;

	if(row >= rows() )
	{
		throw LAException("row >= rows()");
	}
	if(column >= columns() )
	{
		throw LAException("column >= columns()");
	}

	auto predicate = [&](typename Data::Coordinate const & c) { return c.i == row && c.j == column; };

	typename std::list<typename Data::Coordinate>::const_iterator it = std::find_if(data.coordinate_list.begin(), data.coordinate_list.end(), predicate);

	Scalar value = (Scalar) 0.0f;
	if(it != data.coordinate_list.end())
	{
		typename Data::Coordinate c = *it;

		value = c.value;
	}

	return value;
}



template<typename Scalar>
void
SparseCOO<Scalar>::setValue(size_t row, size_t column, Scalar value)
{
	if(row >= rows() )
	{
		throw LAException("SparseCOO::setValue(): row >= rows()" );
	}
	if(column >= columns() )
	{
		throw LAException("SparseCOO::setValue(): column >= columns()");
	}


	this->data.getKeyReference(row,column) = value;
}



template<typename Scalar>
Scalar & 
SparseCOO<Scalar>::operator() (size_t row, size_t column)
{
	if(row >= rows() )
	{
		throw LAException("row >= rows()");
	}
	if(column >= columns() )
	{
		throw LAException("column >= rows()");
	}

	return this->data.getKeyReference(row,column);
}


template<typename Scalar>
void
SparseCOO<Scalar>::resize(size_t rows, size_t columns)
{
	using namespace std;

	data.coordinate_list.clear();

	data.n_rows = rows;
	data.n_columns = columns;
}


template<typename Scalar>
void
SparseCOO<Scalar>::setEye()
{
	data.coordinate_list.clear();
	auto n_elements = rows() < columns()? rows(): columns();

	for(size_t k = 0; k < n_elements; k++)
	{
		typename Data::Coordinate c{k,k,(Scalar)1};
		data.coordinate_list.push_back(c);
	}
}


template<typename Scalar>
size_t 
SparseCOO<Scalar>::nnz() const
{
	return this->data.coordinate_list.size();
}


template<typename Scalar>
void 
SparseCOO<Scalar>::permuteRows(size_t row1, size_t row2)
{
	for(auto it = data.coordinate_list.begin(); it != data.coordinate_list.end(); it++)
	{
		if( it->i == row1)
		{
			it->i = row2;
		}
		else if( it->i == row2)
		{
			it->i = row1;
		}
	}
}


template<typename Scalar>
void 
SparseCOO<Scalar>::permuteColumns(size_t column1, size_t column2)
{
	for(auto it = data.coordinate_list.begin(); it != data.coordinate_list.end(); it++)
	{
		if( it->j == column1)
		{
			it->j = column2;
		}
		else if( it->j == column2)
		{
			it->j = column1;
		}
	}
}



template<typename Scalar>
void 
SparseCOO<Scalar>::symmetricReorder(std::vector<size_t> const & indices)
{
	// check if indices are compatible with the matrix
	if(this->columns() != indices.size())
	{
		throw LAException("indices must match matrix size");
	}

	for(auto it = data.coordinate_list.begin(); it != data.coordinate_list.end(); it++)
	{
		it->i = indices[it->i];
		it->j = indices[it->j];
	}
}



template<typename Scalar>
SparseCOO<Scalar> 
SparseCOO<Scalar>::splitOutBottom(std::size_t row_index)
{
	// validate
	if(row_index >= this->rows())
	{
		throw LAException("row_id is greater than the number of rows");
	}

	SparseCOO<Scalar> bottom_matrix(this->rows()-row_index, this->columns());

	for(typename std::list<typename Data::Coordinate>::iterator it = this->data.coordinate_list.begin(); it != this->data.coordinate_list.end();)
	{
		if( it->i >= row_index)
		{
			size_t new_row_index = it->i-row_index;
			bottom_matrix.setValue(new_row_index, it->j, it->value);
			data.coordinate_list.erase(it++);
		}
		else
		{
			it++;
		}
	}

	// fix row size
	this->data.n_rows = row_index;

	return bottom_matrix;
}


template<typename Scalar>
SparseCOO<Scalar> 
SparseCOO<Scalar>::splitOutRight(std::size_t column_index)
{
	// validate input
	if(column_index >= this->columns())
	{
		throw LAException("column_id is greater than the number of columns");
	}

	SparseCOO<Scalar> right_matrix(this->rows(), this->columns()-column_index);

	for(typename std::list<typename Data::Coordinate>::iterator it = this->data.coordinate_list.begin(); it != this->data.coordinate_list.end();)
	{
		if( it->j >= column_index)
		{
			size_t new_column_index = it->j-column_index;
			right_matrix.setValue(it->i, new_column_index, it->value);
			data.coordinate_list.erase(it++);
		}
		else
		{
			it++;
		}
	}

	// fix column size
	this->data.n_columns = column_index;

	return right_matrix;
}



// Specify traits

template<typename Scalar>
struct Traits<SparseCOO<Scalar> >
{
	static constexpr bool is_writeable() noexcept { return true; }
	static constexpr bool is_resizeable() noexcept { return true; }
	static constexpr bool is_sparse() noexcept { return true; }
};


}	// namespace matrix
}	// namespace mla



#endif
