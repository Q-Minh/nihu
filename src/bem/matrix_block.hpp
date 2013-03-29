/**
 * \file matrix_block.hpp
 * \author Peter Fiala and Peter Rucz, fiala@hit.bme.hu, rucz@hit.bme.hu
 * \brief implementation of a block matrix
 */
 
#ifndef MATRIX_BLOCK_HPP_INCLUDED
#define MATRIX_BLOCK_HPP_INCLUDED

/**
 * \brief proxy class to represent a block of a matrix
 * \details The class is used to represent a matrix block selected by two index vectors.
 * \tparam Matrix the matrix type that is indexed
 * \tparam RowIndex the row index vector type
 * \tparam ColIndex the column index vector type
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
	matrix_block(Matrix &m, RowIndex const &rows, ColIndex const &cols) : m_matrix(m), m_rows(rows), m_cols(cols)
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
 * \brief factory function of matrix_block
 * \details The function allows easy creation of block instances of an arbitrary matrix object
 * \tparam Matrix the matrix type
 * \tparam RowIndex the row index vector type
 * \tparam ColIndex the column index vector type
 * \param matrix [in] the matrix object to index
 * \param rows [in] the rows to be selected
 * \param cols [in] the columns to be selected
 * \return the block proxy
 */
template <class Matrix, class RowIndex, class ColIndex>
matrix_block<Matrix, RowIndex, ColIndex> block(Matrix &matrix, RowIndex const &rows, ColIndex const &cols)
{
	return matrix_block<Matrix, RowIndex, ColIndex>(matrix, rows, cols);
}

#endif // MATRIX_BLOCK_HPP_INCLUDED
