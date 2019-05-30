#include <gtest/gtest.h>

#include "core/field.hpp"
#include "core/function_space.hpp"
#include "core/mesh.hpp"
#include "fmm/elem_center_iterator.hpp"
#include "library/lib_element.hpp"
#include "library/lib_mesh.hpp"
#include "util/casted_iterator.hpp"

template <class Container, class It>
void iterator_test(Container const& elements, It begin, It end)
{
	size_t N = elements.size(), i;

	It first;
	// test preincrement, != and dereference
	for (i = 0, first = begin; first != end; ++first, ++i)
		EXPECT_EQ(*first, elements[i].get_center());

	// test postincrement, == and dereference
	for (i = 0, first = begin; !(first == end); first++, ++i)
		EXPECT_EQ(*first, elements[i].get_center());

	// test indexing
	for (i = 0; i < N; ++i)
		EXPECT_EQ(begin[i], elements[i].get_center());

	// test iterator arithmetics
	EXPECT_EQ(end - begin, N);
	EXPECT_EQ(end - N, begin);
	EXPECT_EQ(begin + N, end);
}

TEST(elem_center_iterator, a)
{
	// create a vector of line elements
	typedef NiHu::line_1_elem elem_t;
	std::vector<elem_t> elements;
	size_t N = 10;
	for (size_t i = 0; i < N; ++i)
	{
		elem_t::coords_t coords;
		coords.setRandom();
		elements.push_back(elem_t(coords));
	}

	iterator_test(elements,
		NiHu::fmm::create_elem_center_iterator(elements.begin()),
		NiHu::fmm::create_elem_center_iterator(elements.end()));

	typedef NiHu::field_view<elem_t, NiHu::field_option::constant> field_t;
	iterator_test(elements,
		NiHu::fmm::create_field_center_iterator(NiHu::casted_iterator<decltype(elements.begin()), field_t>(elements.begin())),
		NiHu::fmm::create_field_center_iterator(NiHu::casted_iterator<decltype(elements.end()), field_t>(elements.end())));
}


TEST(field_center_iterator, a)
{
	using namespace boost::math::double_constants;

	// create circle mesh
	size_t N = 10;
	double R = 1.0;

	auto mesh = NiHu::create_line_1_circle_mesh(R, N);
	auto const& field = NiHu::constant_view(mesh);

	typedef NiHu::line_1_elem elem_t;
	typedef NiHu::field_view<elem_t, NiHu::field_option::constant> field_t;

	auto fb = field.field_begin<field_t>();
	auto fe = field.field_end<field_t>();

	auto fcb = NiHu::fmm::create_field_center_iterator(fb);
	auto fce = NiHu::fmm::create_field_center_iterator(fe);

	std::cout << fce - fcb << std::endl;

	auto elit = mesh.begin<elem_t>();
	for (auto it = fcb; it != fce; ++it)
		EXPECT_EQ(*it, elit++->get_center());

	elit = mesh.begin<elem_t>();
	for (size_t i = 0; i < N; ++i)
		EXPECT_EQ(fcb[i], elit[i].get_center());
}
