#include <mla/matrix/all.h++>
#include <mla/vector/all.h++>


namespace mla {


// explicit template instantiation

// Matrix classes

template class matrix::DenseRowMajor<float>;
template class matrix::DenseRowMajor<double>;

template class matrix::SparseDOK<float>;
template class matrix::SparseDOK<double>;

template class matrix::SparseCOO<float>;
template class matrix::SparseCOO<double>;

template class matrix::Diagonal<float>;
template class matrix::Diagonal<double>;

template class matrix::SparseCRS<float>;
template class matrix::SparseCRS<double>;

template class matrix::SparseCCS<float>;
template class matrix::SparseCCS<double>;

// Vector classes

template class vector::Dense<float>;
template class vector::Dense<double>;

template class vector::SparseCS<float>;
template class vector::SparseCS<double>;


}	// namespace mla

