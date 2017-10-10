#ifndef MLA_MATRIX_TRAITS_HPP
#define MLA_MATRIX_TRAITS_HPP


namespace mla
{
namespace matrix
{

template <class MatrixType>
struct Traits
{
	typedef typename MatrixType::scalar_type scalar_type;

	static constexpr bool is_writeable() noexcept { return false; }
	static constexpr bool is_resizeable() noexcept { return false; }
	static constexpr bool is_sparse() noexcept { return false; }
};


}	// namespace traits
}	// namespace mla

#endif
