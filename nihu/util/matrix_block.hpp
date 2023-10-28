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
 * \file matrix_block.hpp
 * \brief implementation of a block matrix
 * \ingroup util
 */

#ifndef MATRIX_BLOCK_HPP_INCLUDED
#define MATRIX_BLOCK_HPP_INCLUDED

namespace NiHu
{

/**
 * \brief Proxy class to represent a block of a matrix
 * \tparam Matrix the matrix type that is indexed
 * \tparam RowIndex the row index vector type
 * \tparam ColIndex the column index vector type
 * \details
 * The class represents a matrix block selected by two index vectors.
 */
template <class Matrix, class RowIndex, class ColIndex = RowIndex>
class matrix_block
{
public:
	/**
	 * \brief constructor from a matrix and two index vectors
	 * \param [in] m the matrix to be indexed
	 * \param [in] rows the index vector of selected rows
	 * \param [in] cols the index vector of selected columns
	 */
	matrix_block(Matrix &m, RowIndex const &rows, ColIndex const &cols)
		: m_matrix(m), m_rows(rows), m_cols(cols)
	{
	}

	/**
	 * \brief increment the block with a submatrix
	 * \tparam SubMatrix the submatrix type
	 * \param [in] rhs the submatrix to add to the block
	 */
	template <class SubMatrix>
	void operator +=(SubMatrix const &rhs) const
	{
		for (int i = 0; i < m_rows.size(); ++i)
			for (int j = 0; j < m_cols.size(); ++j)
				m_matrix(m_rows(i), m_cols(j)) += rhs(i, j);
	}

protected:
	Matrix &m_matrix;		/**< \brief the matrix to be indexed */
	RowIndex const &m_rows;	/**< \brief the row index vector */
	ColIndex const &m_cols;	/**< \brief the column index vector */
};


/**
 * \brief Factory function of matrix_block
 * \tparam Matrix The matrix type
 * \tparam RowIndex Row index vector type
 * \tparam ColIndex Column index vector type
 * \param[in] matrix The matrix object to index
 * \param[in] rows Rows to be selected
 * \param[in] cols Columns to be selected
 * \return Block proxy to the referenced block
 * 
 * \details 
 * The function allows easy creation of block instances of an arbitrary matrix 
 * object.
 */
template <class Matrix, class RowIndex, class ColIndex>
matrix_block<Matrix, RowIndex, ColIndex>
	block(Matrix &matrix, RowIndex const &rows, ColIndex const &cols)
{
	return matrix_block<Matrix, RowIndex, ColIndex>(matrix, rows, cols);
}

} // end of namespace NiHu

#endif /* MATRIX_BLOCK_HPP_INCLUDED */
