#ifndef MLA_VECTOR_TRAITS_HPP
#define MLA_VECTOR_TRAITS_HPP


namespace mla
{
namespace vector
{

template <class VectorType>
class Traits
{
	typedef typename VectorType::scalar_type scalar_type;

	static constexpr bool is_writeable() noexcept { return false; }
	static constexpr bool is_resizeable() noexcept { return false; }
	static constexpr bool is_sparse() noexcept { return false; }
};


}	// namespace vector
}	// namespace mla

#endif
