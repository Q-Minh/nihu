// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
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
#include <string>
#include <fstream>
#include <stdexcept>

/** \brief Read mesh from OFF format
 * \tparam Tags the element tags to import
 * \param [in] fname the file name
 * \param [in] tags the element tag types
 * \return the imported mesh
 */
template <class...Tags>
mesh<tmp::vector<typename tag2element<Tags>::type...> >
	read_off_mesh(std::string const &fname, Tags...tags)
{
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
	typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

	// read mesh file for reading
	std::ifstream is(fname);
	if (!is)
		throw std::runtime_error("Error reading mesh file");

	// read header from file (first row is 'OFF')
	std::string header;
	if (!(is >> header) || header != "OFF")
		throw std::runtime_error("Possibly invalid off file");

	// read number of nodes and number of elements
	unsigned nNodes, nElements, nEdges;
	if (!(is >> nNodes >> nElements >> nEdges))
		throw std::runtime_error("Error reading number of mesh entries");

	// read nodes
	dMatrix nodes(nNodes, 3);
	for (unsigned i = 0; i < nNodes; ++i)
		if (!(is >> nodes(i,0) >> nodes(i,1) >> nodes(i,2)))
			throw std::runtime_error("Error reading mesh nodes");

	// read elements
	uMatrix elements(nElements, 4);
	for (unsigned i = 0; i < nElements; ++i)
		if (!(is >> elements(i,0) >> elements(i,1) >> elements(i,2) >> elements(i,3)))
			throw std::runtime_error("Error reading mesh elements");

	// close file
	is.close();

	// create and return the mesh
	return create_mesh(nodes, elements, tags...);
}

#endif // READ_OFF_MESH_HPP_INCLUDED
