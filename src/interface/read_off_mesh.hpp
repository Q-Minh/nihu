// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * \file read_off_mesh.hpp
 * \brief export NiHu mesh from OFF format
 */

#ifndef READ_OFF_MESH_HPP_INCLUDED
#define READ_OFF_MESH_HPP_INCLUDED

#include "../core/mesh.hpp"
#include "../library/lib_element.hpp"
#include "../util/type2tag.hpp"

#include <string>
#include <fstream>
#include <stdexcept>

namespace NiHu
{

namespace internal
{
	template <class Tag1, class...Tags>
	struct nvert2elem_id_impl
	{
		typedef typename tag2type<Tag1>::type elem_t;
		static unsigned eval(unsigned nvert)
		{
			if (nvert == elem_t::num_nodes)
				return elem_t::id;
			return nvert2elem_id_impl<Tags...>::eval(nvert);
		}
	};

	template <class Tag1>
	struct nvert2elem_id_impl<Tag1>
	{
		typedef typename tag2type<Tag1>::type elem_t;
		static unsigned eval(unsigned nvert)
		{
			if (nvert == elem_t::num_nodes)
				return elem_t::id;
			throw std::runtime_error("Invalid number of vertices for current element tag vector");
		}
	};
}


/** \brief convert number of vertices to element id
 * \tparam Tags element tags to choose from
 * \param [in] nvert the number of element vertices
 * \param tags the possible element tags
 * \return the first element id from the vector Tags... that matches the number of vertices
 */
template <class...Tags>
unsigned nvert2elem_id(unsigned nvert, Tags...tags)
{
	return internal::nvert2elem_id_impl<Tags...>::eval(nvert);
}


/** \brief Read mesh from OFF format
 * \tparam Tags the element tags to import
 * \param [in] fname the file name
 * \param [in] is the input stream
 * \return the imported mesh
 */
template <class...Tags>
mesh<tmp::vector<typename tag2type<Tags>::type...> >
read_off_mesh(std::istream &is, Tags...tags)
{
	typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

	// read header from file (first row is 'OFF')
	std::string header;
	if (!(is >> header) || header != "OFF")
		throw std::runtime_error("Possibly invalid off file");

	// read number of nodes and number of elements
	unsigned nNodes, nElements, nEdges;
	if (!(is >> nNodes >> nElements >> nEdges))
		throw std::runtime_error("Error reading number of mesh entries");

	// read nodes
	Eigen::MatrixXd nodes(nNodes, 3);
	for (unsigned i = 0; i < nNodes; ++i)
		if (!(is >> nodes(i, 0) >> nodes(i, 1) >> nodes(i, 2)))
			throw std::runtime_error("Error reading mesh nodes");

	// read elements
	uMatrix elements(nElements, 5);
	for (unsigned i = 0; i < nElements; ++i)
	{
		unsigned nvert;
		if (!(is >> nvert))
			throw std::runtime_error("Error reading mesh elements");
		for (unsigned c = 0; c < nvert; ++c)
			if (!(is >> elements(i, c + 1)))
				throw std::runtime_error("Error reading mesh elements");
		elements(i, 0) = nvert2elem_id(nvert, tags...);
	}

	// create and return the mesh
	return create_mesh(nodes, elements, tags...);
}


/** \brief Read mesh from OFF format
 * \tparam Tags the element tags to import
 * \param [in] fname the file name
 * \param [in] tags the element tag types
 * \return the imported mesh
 */
template <class...Tags>
mesh<tmp::vector<typename tag2type<Tags>::type...> >
	read_off_mesh(std::string const &fname, Tags...tags)
{
	typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

	// read mesh file for reading
	std::ifstream is(fname);
	if (!is)
		throw std::runtime_error("Error reading mesh file");

	return read_off_mesh(is, tags...);
}

}

#endif // READ_OFF_MESH_HPP_INCLUDED

