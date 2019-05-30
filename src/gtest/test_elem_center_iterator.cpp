#include <gtest/gtest.h>

#include <boost/math/constants/constants.hpp>

#include "core/field.hpp"
#include "core/function_space.hpp"
#include "core/mesh.hpp"
#include "fmm/elem_center_iterator.hpp"
#include "library/lib_element.hpp"
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

	Eigen::Matrix<double, 2, Eigen::Dynamic> nodes(2, N);
	for (size_t i = 0; i < N; ++i)
	{
		double phi = two_pi / N * i;
		nodes(0, i) = R * std::cos(phi);
		nodes(1, i) = R * std::sin(phi);
	}

	Eigen::Matrix<unsigned, 2 + 1, Eigen::Dynamic> elements(3, N);
	for (size_t i = 0; i < N; ++i)
	{
		elements(0, i) = NiHu::line_1_elem::id;
		elements(1, i) = i;
		elements(2, i) = (i + 1) % N;
	}

	typedef tmp::vector<NiHu::line_1_elem> elem_type_vector_t;
	typedef NiHu::mesh<elem_type_vector_t> mesh_t;
	typedef NiHu::field_view<NiHu::line_1_elem, NiHu::field_option::constant> field_t;
	typedef NiHu::function_space_view<mesh_t, NiHu::field_option::constant> function_space_t;

	mesh_t mesh;
	for (size_t i = 0; i < N; ++i)
		mesh.add_node(nodes.col(i).data());

	for (size_t i = 0; i < N; ++i)
	{
		unsigned el[3];
		el[0] = elements(0, i);
		el[1] = elements(1, i);
		el[2] = elements(2, i);
		mesh.add_elem(el);
	}
	function_space_t const& field = NiHu::constant_view(mesh);

	auto fb = field.field_begin<field_t>();
	auto fe = field.field_end<field_t>();

	auto fcb = NiHu::fmm::create_field_center_iterator(fb);
	auto fce = NiHu::fmm::create_field_center_iterator(fe);

	std::cout << fce - fcb << std::endl;

	for (auto it = fcb; it != fce; ++it)
		EXPECT_EQ(*it, (nodes.col(it - fcb) + nodes.col((it - fcb + 1) % N)) / 2);

	for (size_t i = 0; i < N; ++i)
		EXPECT_EQ(fcb[i], (nodes.col(i) + nodes.col((i + 1) % N)) / 2);
}
