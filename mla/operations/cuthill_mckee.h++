#ifndef MLA_OPERATIONS_CUTHILL_MCKEE_HPP
#define MLA_OPERATIONS_CUTHILL_MCKEE_HPP

#include <type_traits>
#include <list>
#include <vector>
#include <algorithm>

#include <mla/LAException.h++>

#include <mla/matrix/all.h++>


namespace mla {

/**
 * Implementation of Cuthill-McKee algorithms
 *
 **/

template<typename Scalar>
std::vector<size_t>
cuthill_mckee(mla::matrix::SparseDOK<Scalar> const &A)
{
	// matrix A must be a square matrix
	if( !A.isSquare() )
	{
		throw LAException("A must be a square matrix");
	}
	
	size_t A_size = A.rows();

	//calculate the degree of each vertex
	class Element {
		size_t m_index;
		size_t m_degree;

	public:
		Element(size_t index, size_t degree = 0) 
			: m_index(index), m_degree(degree)
		{ };

		/**
		 * order operator that imposes the Cuthill-McKee order
		 */
		bool operator<(Element const &rhs) const
		{
			return this->m_degree > rhs.getDegree();
		}

		bool operator==(Element const &rhs) const
		{
			if (this->m_index != rhs.getIndex()) 
				return false;
			if (this->m_degree != rhs.getDegree())
				return false;
			return true;
		}

		/**
		 * helper function that increments the degree
		 **/
		void incrementDegree() {this->m_degree++;};

		size_t getIndex() const {return this->m_index;};
		size_t getDegree() const {return this->m_degree;};
	};


	std::list<Element> adjacency;

	// init adjacency data structure
	for(size_t i = 0; i < A_size; i++)
	{
		adjacency.push_back( Element(i) );
	}

	// calculate the degree of each vertex
	for(auto &kv: A.data.key_value_map)
	{
		size_t key = kv.first.first;	// kv.first is object of type pair<size_t,size_t>(i,j)

		auto compare_predicate = [&](Element const &e) {return e.getIndex() == key;};
		typename std::list<Element>::iterator it = std::find_if(adjacency.begin(), adjacency.end(), compare_predicate);

		if(it == adjacency.end())
		{
			throw LAException("Key doesn''t exist");
		}

		it->incrementDegree();
	}
	
	// reorder adjacency list based on degree
	adjacency.sort();

	// assign index to vector to output
	std::vector<size_t> index;
	index.reserve(adjacency.size());
	for(auto &e: adjacency)
	{
		index.push_back( e.getIndex() );
	}


	return index;
}


}	// mla

#endif
