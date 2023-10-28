/// \file m2l_indices.hpp
/// \brief implementation of class template fmm::m2l_indices

#ifndef M2L_INDICES_HPP_INCLUDED
#define M2L_INDICES_HPP_INCLUDED

#include "bounding_box.hpp"
#include <type_traits>

namespace NiHu
{
namespace fmm
{
/// \brief class assigning indices to M2L distances
/// \tparam Dim the space dimensionality
template <unsigned Dim>
class m2l_indices
{
public:
	static size_t const dimension = Dim;
	typedef bounding_box<dimension> bounding_box_t;
	typedef typename bounding_box_t::location_t location_t;

private:
	static size_t drel2idx(location_t const &drel, std::integral_constant<size_t, 1U>)
	{
		return size_t(round(drel(0)) + 3);
	}

	static size_t drel2idx(location_t const &drel, std::integral_constant<size_t, 2U>)
	{
		return size_t(round(drel(0)) + 3 + 7 * (round(drel(1)) + 3));
	}

	static size_t drel2idx(location_t const &drel, std::integral_constant<size_t, 3U>)
	{
		return size_t(round(drel(0)) + 3 + 7 * (round(drel(1)) + 3) + 7 * 7 * (round(drel(2)) + 3));
	}
public:

	static size_t eval(bounding_box_t const &to, bounding_box_t const &from)
	{
		location_t drel = (to.get_center() - from.get_center()) / to.get_diameter();
		return drel2idx(drel, std::integral_constant<size_t, dimension>());
	}
};
} // end of namespace fmm
} // namespace NiHu

#endif
