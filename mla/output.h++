#ifndef MLA_OUTPUT_HPP
#define MLA_OUTPUT_HPP

#include <iostream>


namespace mla
{


#ifdef MLA_MATRIX_HPP
template<typename scalar, template<typename> class MatrixStoragePolicy>
std::ostream & operator << ( std::ostream &out, Matrix<scalar, MatrixStoragePolicy> &m)
{
	using namespace std;

	for(size_t i = 0; i < m.rows(); i++)
	{
		out << "[\t" << m.getValue(i,0);
		for(size_t j = 1; j < m.columns(); j++)
		{
			out << ",\t" << m.getValue(i,j);
		}
		out << "]" << endl;
	}
	return out;
}


template<typename scalar, template<typename> class MatrixStoragePolicy>
void dump_octave( std::ostream &out, std::string name, Matrix<scalar, MatrixStoragePolicy> &m)
{
	using namespace std;

	out << "# Created by mla\n";
	out << "# name: " << name << "\n";
	out << "# type: matrix\n";
	out << "# rows: " << m.rows() << "\n";
	out << "# columns: " << m.columns() << "\n";

	for(size_t j = 0; j < m.columns(); j++)
	{
		for(size_t i = 0; i < m.rows(); i++)
		{
			out << " " << m.getValue(i,j);
		}
		out << "\n";
	}
	out << endl;
}


template<typename scalar, template<typename> class VectorStoragePolicy>
void dump_octave( std::ostream &out, std::string name, Vector<scalar, VectorStoragePolicy> &v)
{
	using namespace std;

	out << "# Created by mla\n";
	out << "# name: " << name << "\n";
	out << "# type: matrix\n";
	out << "# rows: " << v.size() << "\n";
	out << "# columns: 1\n";
	for(size_t i = 0; i < v.size(); i++)
	{
		out << " " << v.getValue(i) << "\n";
	}
	out << endl;
}

#endif


#ifdef MLA_VECTOR_HPP
template<typename scalar, template<typename> class VectorStoragePolicy>
std::ostream & operator << ( std::ostream &out, Vector<scalar, VectorStoragePolicy> v)
{
	using namespace std;

	for(size_t i = 0; i < v.size(); i++)
	{
		out << "[\t" << v.getValue(i) << "\t]" << endl;
	}
	return out;
}


#endif


}	// namespace mla

#endif
